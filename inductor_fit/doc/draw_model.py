# Draws the inductor_fit model topologies (doc/inductor_fit_model.png for 2-port inductors,
# doc/inductor_fit_ct_model.png for center-tapped inductors) with schemdraw
# (pip install schemdraw). Only needed to regenerate the documentation images.

import os
import schemdraw
import schemdraw.elements as elm
from schemdraw.util import Point

U = 2.4  # element length
DOC_DIR = os.path.dirname(os.path.abspath(__file__))

def line(d, a, b, **kwargs):
    d.add(elm.Line(**kwargs).endpoints(a, b))

def parallel_horizontal(d, a, top, bottom, gap=1.0):
    # 'top' and 'bottom' elements in parallel from point a to the right, returns the right node
    a = Point(a)
    b = a + Point((U + 1.0, 0))
    line(d, a, a + Point((0.5, 0)))
    line(d, b - Point((0.5, 0)), b)
    for elem, dy in ((top, gap / 2), (bottom, -gap / 2)):
        p, q = a + Point((0.5, dy)), b + Point((-0.5, dy))
        line(d, a + Point((0.5, 0)), p)
        d.add(elem.endpoints(p, q))
        line(d, q, b + Point((-0.5, 0)))
    return b

def shunt_network(d, top, cox_label, rsi_label, csi_label, node_label, side):
    # Cox down from 'top' to the substrate node, then Rsi || Csi to ground; returns the substrate node
    s = Point(top) - Point((0, U + 1.0))
    loc = lambda side: 'top' if side == 'left' else 'bottom'  # label sides of downward elements
    d.add(elm.Capacitor().endpoints(top, s).label(cox_label, loc(side)))
    d.add(elm.Dot().at(s).label(node_label, 'right' if side == 'left' else 'left', ofst=0.2))
    w = 0.8
    for elem, dx, label_side in ((elm.Resistor(), -w, 'left'), (elm.Capacitor(), w, 'right')):
        p = s + Point((dx, -0.6))
        q = p - Point((0, U))
        line(d, s, s + Point((dx, 0)))
        line(d, s + Point((dx, 0)), p)
        d.add(elem.endpoints(p, q).label(rsi_label if dx < 0 else csi_label, loc(label_side)))
        line(d, q, q - Point((0, 0.4)))
        line(d, q - Point((0, 0.4)), s + Point((0, -U - 1.0)))
    d.add(elm.Ground().at(s + Point((0, -U - 1.0))))
    return s

def series_coil(d, start, suffix=''):
    # Rs, Ls and two skin sections from 'start' to the right; returns the end point and the
    # center of Ls (for the coupling bracket)
    a = start + Point((U, 0))
    d.add(elm.Resistor().endpoints(start, a).label('Rs' + suffix))
    b = a + Point((U, 0))
    d.add(elm.Inductor2(loops=3).endpoints(a, b).label('Ls' + suffix))
    c = parallel_horizontal(d, b, elm.Inductor2(loops=2).label('Lskin1' + suffix), elm.Resistor().label('Rskin1' + suffix, 'bottom'))
    e = parallel_horizontal(d, c, elm.Inductor2(loops=2).label('Lskin2' + suffix), elm.Resistor().label('Rskin2' + suffix, 'bottom'))
    return e, Point(((a.x + b.x) / 2, a.y))

def cs_across(d, n1, n2, h):
    line(d, n1, n1 + Point((0, h)))
    line(d, n2, n2 + Point((0, h)))
    mid = (n1.x + n2.x) / 2
    line(d, n1 + Point((0, h)), Point((mid - U / 2, n1.y + h)))
    d.add(elm.Capacitor().endpoints(Point((mid - U / 2, n1.y + h)), Point((mid + U / 2, n1.y + h))).label('Cs'))
    line(d, Point((mid + U / 2, n1.y + h)), n2 + Point((0, h)))
    return mid

def save(d, filename):
    d.save(os.path.join(DOC_DIR, filename), dpi=110, transparent=False)


def draw_two_port():
    d = schemdraw.Drawing(show=False)
    d.config(unit=U, fontsize=12)

    # ---- series branch ----
    p1 = Point((0, 0))
    d.add(elm.Dot(open=True).at(p1).label('p1', 'left'))
    n1 = p1 + Point((1.5, 0))
    line(d, p1, n1)
    d.add(elm.Dot().at(n1))
    e, _ = series_coil(d, n1)
    n2 = e + Point((0.5, 0))
    line(d, e, n2)
    d.add(elm.Dot().at(n2))
    p2 = n2 + Point((1.5, 0))
    line(d, n2, p2)
    d.add(elm.Dot(open=True).at(p2).label('p2', 'right'))

    # ---- Cs across the coil, substrate networks ----
    mid = cs_across(d, n1, n2, 2.0)
    s1 = shunt_network(d, n1, 'Cox1', 'Rsi1', 'Csi1', 's1', 'left')
    s2 = shunt_network(d, n2, 'Cox2', 'Rsi2', 'Csi2', 's2', 'right')

    # ---- substrate coupling s1 - s2 (only used if needed) ----
    w = U + 1.0
    left = Point((mid - w / 2, s1.y))
    line(d, s1 + Point((0.8, 0)), left)
    right = parallel_horizontal(d, left, elm.Resistor().label('Rsub12'), elm.Capacitor().label('Csub12', 'bottom'), gap=1.4)
    line(d, right, s2 - Point((0.8, 0)))

    save(d, 'inductor_fit_model.png')


def draw_center_tap():
    d = schemdraw.Drawing(show=False)
    d.config(unit=U, fontsize=12)

    # ---- half coils p1 -> m -> p2 ----
    p1 = Point((0, 0))
    d.add(elm.Dot(open=True).at(p1).label('p1', 'left'))
    n1 = p1 + Point((1.5, 0))
    line(d, p1, n1)
    d.add(elm.Dot().at(n1))
    e1, ls1 = series_coil(d, n1, '_h1')
    m = e1 + Point((0.5, 0))
    line(d, e1, m)
    d.add(elm.Dot().at(m).label('m', 'top', ofst=0.15))
    m2 = m + Point((2.2, 0))
    line(d, m, m2)
    d.add(elm.Dot().at(m2))
    e2, ls2 = series_coil(d, m2 + Point((0.5, 0)), '_h2')
    line(d, m2, m2 + Point((0.5, 0)))
    n2 = e2 + Point((0.5, 0))
    line(d, e2, n2)
    d.add(elm.Dot().at(n2))
    p2 = n2 + Point((1.5, 0))
    line(d, n2, p2)
    d.add(elm.Dot(open=True).at(p2).label('p2', 'right'))

    # ---- magnetic coupling of Ls_h1 and Ls_h2 ----
    hk = 1.6
    for ls in (ls1, ls2):
        line(d, ls + Point((0, 1.0)), ls + Point((0, hk)), ls='--')
    line(d, ls1 + Point((0, hk)), ls2 + Point((0, hk)), ls='--')
    d.add(elm.Label().at(Point(((ls1.x + ls2.x) / 2, hk + 0.35))).label('k (Ls_h1 coupled to Ls_h2)'))

    # ---- Cs across the coil ----
    cs_across(d, n1, n2, 2.8)

    # ---- center tap lead ----
    ct = m - Point((0, U + 1.0))
    d.add(elm.Resistor().endpoints(m, ct).label('Rct', 'top'))
    d.add(elm.Dot(open=True).at(ct).label('ct', 'bottom'))

    # ---- substrate networks at p1, m, p2 ----
    shunt_network(d, n1, 'Cox1', 'Rsi1', 'Csi1', 's1', 'left')
    shunt_network(d, m2, 'Coxct', 'Rsict', 'Csict', 'sct', 'right')
    shunt_network(d, n2, 'Cox2', 'Rsi2', 'Csi2', 's2', 'right')

    save(d, 'inductor_fit_ct_model.png')


draw_two_port()
draw_center_tap()
