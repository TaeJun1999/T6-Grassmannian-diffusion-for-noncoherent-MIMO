#!/usr/bin/env python3
"""T6 scooping check via arXiv OAI-PMH (bulk metadata incl. abstracts; separate rate limit
from the search API).  Harvests every record in the given sets modified since --since,
keeps those whose categories intersect --cats (or whose text contains a core topic word),
and scores title+abstract with tiered keyword lists.

Usage:
    python3 t6_arxiv_oai_check.py --since 2026-08-25 --sets cs,eess \
        --cats cs.IT,eess.SP,cs.LG --out t6_oai_2026-09-16.json
"""
import argparse, json, re, time, datetime, sys
import urllib.request, urllib.parse, urllib.error, xml.etree.ElementTree as ET

UA = "t6-scooping-check/1.0"
OAI = "https://oaipmh.arxiv.org/oai"
NS = {'o': 'http://www.openarchives.org/OAI/2.0/', 'a': 'http://arxiv.org/OAI/arXiv/'}

CORE = ['grassmann', 'stiefel', 'noncoherent', 'non-coherent', 'unitary space-time',
        'unitary space time', 'partially coherent', 'bingham', 'ustm']
TOPIC = CORE + ['constellation', 'subspace', 'blind detection', 'pilot-free', 'pilotless',
                'without csi', 'no csi', 'csi-free', 'differential detection', 'noncoherent',
                'block fading', 'block-fading']
GEN = ['diffusion model', 'diffusion prior', 'diffusion-based', 'score-based', 'score based',
       'generative', 'flow matching', 'flow-matching', 'denoising', 'learned prior',
       'learnable prior', 'deep prior', 'generative prior']
METHOD = GEN + ['riemannian', 'manifold', 'prior', 'learned', 'learning', 'neural', 'deep']
COMM = ['mimo', 'wireless', 'channel', 'antenna', 'fading', 'detection', 'receiver', 'modulation']


def fetch(url, tries=5):
    req = urllib.request.Request(url, headers={'User-Agent': UA})
    for k in range(tries):
        try:
            return urllib.request.urlopen(req, timeout=180).read()
        except urllib.error.HTTPError as ex:
            wait = int(ex.headers.get('Retry-After', '20')) if ex.code in (503, 429) else 20
            if k == tries - 1:
                raise
            print(f"   HTTP {ex.code}, waiting {wait}s", flush=True)
            time.sleep(wait)
        except Exception as ex:
            if k == tries - 1:
                raise
            print(f"   {ex!r}, retry", flush=True)
            time.sleep(15)


def parse(xml_bytes):
    root = ET.fromstring(xml_bytes)
    recs = []
    for r in root.iter('{%s}record' % NS['o']):
        m = r.find('.//a:arXiv', NS)
        if m is None:  # deleted record
            continue
        g = lambda tag: (m.findtext('a:' + tag, default='', namespaces=NS) or '').strip()
        authors = ['%s %s' % ((a.findtext('a:forenames', default='', namespaces=NS) or '').strip(),
                              (a.findtext('a:keyname', default='', namespaces=NS) or '').strip())
                   for a in m.findall('a:authors/a:author', NS)]
        recs.append(dict(id=g('id'), created=g('created'), updated=g('updated'),
                         cats=g('categories').split(), title=' '.join(g('title').split()),
                         abstract=' '.join(g('abstract').split()), authors=authors,
                         datestamp=(r.findtext('o:header/o:datestamp', default='', namespaces=NS) or '')))
    tok = root.find('.//o:resumptionToken', NS)
    return recs, (tok.text.strip() if tok is not None and tok.text else None)


def harvest(set_, since, sleep):
    url = OAI + '?' + urllib.parse.urlencode({'verb': 'ListRecords', 'metadataPrefix': 'arXiv',
                                              'set': set_, 'from': since})
    recs, page = [], 0
    while True:
        page += 1
        rs, tok = parse(fetch(url))
        recs += rs
        print(f"   set={set_} page={page} records={len(rs)} total={len(recs)}", flush=True)
        if not tok:
            break
        url = OAI + '?' + urllib.parse.urlencode({'verb': 'ListRecords', 'resumptionToken': tok})
        time.sleep(sleep)
    return recs


def hits(text, words):
    t = text.lower()
    # word-boundary match (no letters immediately before/after) so that e.g. 'ustm'
    # does not match 'adjustment' and 'prior' does not match 'priority'
    # trailing boundary only for short tokens (<= 5 chars); stems like 'grassmann' must
    # still match 'grassmannian(s)'
    return sorted({w for w in words if re.search(
        r'(?<![a-z])' + re.escape(w) + (r'(?![a-z])' if len(w) <= 5 else ''), t)})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--since', default='2026-08-25')
    ap.add_argument('--sets', default='cs,eess')
    ap.add_argument('--cats', default='cs.IT,eess.SP,cs.LG,math.IT')
    ap.add_argument('--sleep', type=float, default=3.0)
    ap.add_argument('--out', default='t6_oai_check.json')
    args = ap.parse_args()
    cats = set(args.cats.split(','))

    pool = {}
    for s in args.sets.split(','):
        for r in harvest(s, args.since, args.sleep):
            pool[r['id']] = r
        time.sleep(args.sleep)
    print(f"harvested unique records: {len(pool)}", flush=True)
    import gzip
    with gzip.open(args.out.replace('.json', '') + '_pool.json.gz', 'wt') as f:
        json.dump(list(pool.values()), f, ensure_ascii=False)

    rows = []
    for r in pool.values():
        text = r['title'] + ' ' + r['abstract']
        core, topic, gen, meth, comm = (hits(text, CORE), hits(text, TOPIC), hits(text, GEN),
                                        hits(text, METHOD), hits(text, COMM))
        in_cats = bool(set(r['cats']) & cats)
        if not (in_cats or core):
            continue
        if core and gen and comm:
            level = 'RED'       # core topic + generative model + communications
        elif core and gen:
            level = 'ORANGE'    # core topic + generative model (any field)
        elif core and comm and meth:
            level = 'YELLOW'    # core topic + comms + some learning/manifold word
        elif topic and gen and comm:
            level = 'BLUE'      # weak topic word + generative + comms
        elif core and comm:
            level = 'GREEN'     # core topic + comms, classical (baselines to watch)
        else:
            continue
        rows.append(dict(id=r['id'], created=r['created'], updated=r['updated'], datestamp=r['datestamp'],
                         cats=r['cats'], title=r['title'], authors=r['authors'][:8],
                         abstract=r['abstract'], level=level, core=core, topic=topic, gen=gen,
                         method=meth, comm=comm, new=(r['created'] >= args.since)))
    order = {'RED': 0, 'ORANGE': 1, 'YELLOW': 2, 'BLUE': 3, 'GREEN': 4}
    rows.sort(key=lambda p: (order[p['level']], p['id']))
    n_in_cats = sum(1 for r in pool.values() if set(r['cats']) & cats)

    out = dict(run_date=datetime.date.today().isoformat(), since=args.since, sets=args.sets.split(','),
               cats=sorted(cats), n_harvested=len(pool), n_in_cats=n_in_cats,
               tiers=dict(CORE=CORE, TOPIC=TOPIC, GEN=GEN, METHOD=METHOD, COMM=COMM),
               counts={k: sum(1 for p in rows if p['level'] == k) for k in order},
               candidates=rows)
    with open(args.out, 'w') as f:
        json.dump(out, f, indent=1, ensure_ascii=False)

    print(f"\nrecords in target categories: {n_in_cats}; levels: {out['counts']}")
    print("\n==== CANDIDATES ====")
    for p in rows:
        if p['level'] == 'GREEN' and len(rows) > 60:
            continue
        print(f"\n[{p['level']}] arXiv:{p['id']} {'NEW' if p['new'] else 'UPD'} created={p['created']} "
              f"updated={p['updated']} {p['cats']}")
        print(f"   {p['title']}")
        print(f"   {', '.join(p['authors'][:6])}")
        print(f"   core={p['core']} gen={p['gen']} comm={p['comm']}")
        print(f"   {p['abstract'][:400]}")
    print("\nDONE", flush=True)


if __name__ == '__main__':
    main()
