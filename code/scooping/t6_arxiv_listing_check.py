#!/usr/bin/env python3
"""T6 scooping check via arXiv monthly listing pages (fallback when the API returns 429).

Stage 1: download https://arxiv.org/list/<cat>/<YYYY-MM>?skip=..&show=2000 for the
         categories/months given, parse (id, title, authors, subjects), dedup.
Stage 2: title-level keyword scan (tier A = topic words, tier B = method words).
Stage 3: fetch https://arxiv.org/abs/<id> for tier-A title matches (plus tier-B
         matches in cs.IT/eess.SP) and scan the abstract for the other tier.

Usage:
    python3 t6_arxiv_listing_check.py --months 2026-08,2026-09 \
        --full_cats cs.IT,eess.SP --title_cats cs.LG --out t6_listing_2026-09-16.json
"""
import argparse, json, re, time, datetime, html
import urllib.request

UA = "Mozilla/5.0 (compatible; t6-scooping-check/1.0)"

TIER_A = ['grassmann', 'stiefel', 'noncoherent', 'non-coherent', 'unitary space-time',
          'unitary space time', 'constellation', 'subspace', 'partially coherent',
          'bingham', 'blind detection', 'pilot-free', 'pilotless', 'without csi', 'no csi']
TIER_B = ['diffusion', 'generative', 'flow matching', 'flow-matching', 'score-based',
          'score based', 'riemannian', 'manifold', 'learned prior', 'learnable prior',
          'deep prior', 'prior']
TIER_C = ['learn', 'neural', 'deep', 'transformer', 'autoencoder']  # weak method words


def fetch(url, tries=3, wait=10):
    req = urllib.request.Request(url, headers={'User-Agent': UA})
    for k in range(tries):
        try:
            return urllib.request.urlopen(req, timeout=120).read().decode('utf-8', 'replace')
        except Exception as ex:
            if k == tries - 1:
                raise
            print(f"   retry {k+1} on {url}: {ex!r}", flush=True)
            time.sleep(wait)


ENTRY_RE = re.compile(
    r'<dt>.*?<a href\s*=\s*"/abs/(?P<id>[^"]+)"\s+title="Abstract"[^>]*>.*?</dt>\s*<dd>(?P<body>.*?)</dd>',
    re.S)
TITLE_RE = re.compile(r'<div class=[\'"]list-title[^\'"]*[\'"]>\s*<span class=[\'"]descriptor[\'"]>Title:</span>(.*?)</div>', re.S)
AUTH_RE = re.compile(r'<div class=[\'"]list-authors[\'"]>(.*?)</div>', re.S)
SUBJ_RE = re.compile(r'<div class=[\'"]list-subjects[\'"]>\s*<span class=[\'"]descriptor[\'"]>Subjects:</span>(.*?)</div>', re.S)
TAG_RE = re.compile(r'<[^>]+>')
TOTAL_RE = re.compile(r'Total of (\d+) entries')


def clean(s):
    return ' '.join(html.unescape(TAG_RE.sub('', s)).split())


def parse_listing(page):
    out = []
    for m in ENTRY_RE.finditer(page):
        body = m.group('body')
        t, a, s = TITLE_RE.search(body), AUTH_RE.search(body), SUBJ_RE.search(body)
        out.append(dict(id=m.group('id'),
                        title=clean(t.group(1)) if t else '',
                        authors=clean(a.group(1)) if a else '',
                        subjects=clean(s.group(1)) if s else ''))
    tot = TOTAL_RE.search(page)
    return out, (int(tot.group(1)) if tot else None)


def get_month(cat, ym, sleep):
    entries, skip, total = [], 0, None
    while True:
        url = f"https://arxiv.org/list/{cat}/{ym}?skip={skip}&show=2000"
        page = fetch(url)
        es, total = parse_listing(page)
        entries += es
        print(f"   {cat} {ym}: skip={skip} got={len(es)} total={total}", flush=True)
        skip += 2000
        if not es or total is None or skip >= total:
            break
        time.sleep(sleep)
    return entries, total


ABS_RE = re.compile(r'<blockquote class="abstract[^"]*">\s*<span class="descriptor">Abstract:</span>(.*?)</blockquote>', re.S)
DATE_RE = re.compile(r'<div class="dateline">(.*?)</div>', re.S)


def get_abstract(aid):
    page = fetch(f"https://arxiv.org/abs/{aid}")
    m, d = ABS_RE.search(page), DATE_RE.search(page)
    return (clean(m.group(1)) if m else ''), (clean(d.group(1)) if d else '')


def hits(text, words):
    t = text.lower()
    return [w for w in words if w in t]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--months', default='2026-08,2026-09')
    ap.add_argument('--full_cats', default='cs.IT,eess.SP')
    ap.add_argument('--title_cats', default='cs.LG')
    ap.add_argument('--sleep', type=float, default=2.0)
    ap.add_argument('--out', default='t6_listing_check.json')
    ap.add_argument('--skip_abstracts', action='store_true')
    args = ap.parse_args()
    months = args.months.split(',')
    full_cats = [c for c in args.full_cats.split(',') if c]
    title_cats = [c for c in args.title_cats.split(',') if c]

    pool, listing_stats = {}, []
    for cat in full_cats + title_cats:
        for ym in months:
            es, tot = get_month(cat, ym, args.sleep)
            listing_stats.append(dict(cat=cat, month=ym, total=tot, parsed=len(es)))
            for e in es:
                p = pool.setdefault(e['id'], dict(e, cats=[]))
                p['cats'].append(cat)
            time.sleep(args.sleep)
    print(f"pooled unique entries: {len(pool)}", flush=True)

    # Stage 2: title scan
    for p in pool.values():
        p['title_A'] = hits(p['title'], TIER_A)
        p['title_B'] = hits(p['title'], TIER_B)
        p['title_C'] = hits(p['title'], TIER_C)
    A = [p for p in pool.values() if p['title_A']]
    AB = [p for p in A if p['title_B']]
    print(f"title tier-A matches: {len(A)}; A and B: {len(AB)}", flush=True)

    # Stage 3: abstracts for: all tier-A titles; tier-B titles in the full cats;
    # tier-A words 'grassmann'/'stiefel' anywhere (already in A).
    cand = {p['id']: p for p in A}
    for p in pool.values():
        if p['title_B'] and any(c in full_cats for c in p['cats']):
            cand[p['id']] = p
    print(f"fetching abstracts for {len(cand)} candidates", flush=True)
    if not args.skip_abstracts:
        for i, (aid, p) in enumerate(sorted(cand.items())):
            try:
                ab, dl = get_abstract(aid)
            except Exception as ex:
                ab, dl = '', f'ERROR {ex!r}'
            p['abstract'], p['dateline'] = ab, dl
            p['abs_A'], p['abs_B'], p['abs_C'] = hits(ab, TIER_A), hits(ab, TIER_B), hits(ab, TIER_C)
            if (i + 1) % 10 == 0:
                print(f"   abstracts {i+1}/{len(cand)}", flush=True)
            time.sleep(args.sleep)

    # Scoring: strong = (A in title or abstract) and (B in title or abstract)
    rows = []
    for p in cand.values():
        a = set(p['title_A']) | set(p.get('abs_A', []))
        b = set(p['title_B']) | set(p.get('abs_B', []))
        c = set(p['title_C']) | set(p.get('abs_C', []))
        core = {'grassmann', 'stiefel', 'noncoherent', 'non-coherent', 'unitary space-time',
                'unitary space time', 'partially coherent', 'bingham'} & a
        strongB = ({'diffusion', 'generative', 'flow matching', 'flow-matching', 'score-based',
                    'score based', 'learned prior', 'learnable prior', 'deep prior'} & b)
        if core and strongB:
            level = 'RED'      # topic core + generative-model word
        elif core and (b or c):
            level = 'ORANGE'   # topic core + any method word
        elif core:
            level = 'YELLOW'   # topic core only (baseline / constellation design)
        elif a and strongB:
            level = 'BLUE'     # weak topic word + generative-model word
        else:
            level = 'GREY'
        p['level'] = level
        rows.append(p)
    order = {'RED': 0, 'ORANGE': 1, 'YELLOW': 2, 'BLUE': 3, 'GREY': 4}
    rows.sort(key=lambda p: (order[p['level']], p['id']))

    out = dict(run_date=datetime.date.today().isoformat(), months=months, full_cats=full_cats,
               title_cats=title_cats, listing_stats=listing_stats, n_pool=len(pool),
               tiers=dict(A=TIER_A, B=TIER_B, C=TIER_C), candidates=rows)
    with open(args.out, 'w') as f:
        json.dump(out, f, indent=1, ensure_ascii=False)

    print("\n==== CANDIDATES (RED/ORANGE/YELLOW/BLUE) ====")
    for p in rows:
        if p['level'] == 'GREY':
            continue
        print(f"\n[{p['level']}] arXiv:{p['id']}  {p['cats']}  {p.get('dateline','')}")
        print(f"   {p['title']}")
        print(f"   {p['authors'][:160]}")
        print(f"   subjects: {p['subjects'][:120]}")
        print(f"   A={sorted(set(p['title_A'])|set(p.get('abs_A',[])))} B={sorted(set(p['title_B'])|set(p.get('abs_B',[])))}")
        if p.get('abstract'):
            print(f"   abstract: {p['abstract'][:500]}")
    print("\nDONE", flush=True)


if __name__ == '__main__':
    main()
