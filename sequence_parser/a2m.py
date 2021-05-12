from itertools import groupby

def read_a2m(file):
    """
    generator: yields entries in a a2m file, (identifier, sequence-feature list)

    ignores periods
    
    """
    is_head = lambda x: x.startswith('>')
    parse_h = lambda x: tuple(x)[-1].strip()[1:]
    parse_s = lambda x: ''.join(j for i in x for j in i if j!='.' and not j.isspace())
    it      = iter(groupby(open(file), is_head))
    it      = iter(groupby(open(file), is_head)) if next(it)[0] else it
    for k, g in it:
        if k:
            h = parse_h(g)
        else:
            yield h, parse_s(g)
