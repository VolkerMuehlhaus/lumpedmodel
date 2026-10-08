# Draws the inductor_fit model topology (doc/inductor_fit_model.png) with schemdraw
# (pip install schemdraw). Only needed to regenerate the documentation image.

import os
import schemdraw
import schemdraw.elements as elm
from schemdraw.util import Point

U = 2.4  # element length

d = schemdraw.Drawing(show=False)
d.config(unit=U, fontsize=12)

def line(a, b):
    d.add(elm.Line().endpoints(a, b))

def parallel_horizontal(a, top, bottom, gap=1.0):
    # 'top' and 'bottom' elements in parallel from point a to the right, returns the right node
    a = Point(a)
    b = a + Point((U + 1.0, 0))
    line(a, a + Point((0.5, 0)))
    line(b - Point((0.5, 0)), b)
    for elem, dy in ((top, gap / 2), (bottom, -gap / 2)):
        p, q = a + Point((0.5, dy)), b + Point((-0.5, dy))
        line(a + Point((0.5, 0)), p)
        d.add(elem.endpoints(p, q))
        line(q, b + Point((-0.5, 0)))
    return b

def shunt_network(top, cox_label, rsi_label, csi_label, node_label, side):
    # Cox down from 'top' to the substrate node, then Rsi || Csi to ground; returns the substrate node
    s = Point(top) - Point((0, U + 1.0))
    loc = lambda side: 'top' if side == 'left' else 'bottom'  # label sides of downward elements
    d.add(elm.Capacitor().endpoints(top, s).label(cox_label, loc(side)))
    d.add(elm.Dot().at(s).label(node_label, 'right' if side == 'left' else 'left', ofst=0.2))
    w = 0.8
    for elem, dx, label_side in ((elm.Resistor(), -w, 'left'), (elm.Capacitor(), w, 'right')):
        p = s + Point((dx, -0.6))
        q = p - Point((0, U))
        line(s, s + Point((dx, 0)))
        line(s + Point((dx, 0)), p)
        d.add(elem.endpoints(p, q).label(rsi_label if dx < 0 else csi_label, loc(label_side)))
        line(q, q - Point((0, 0.4)))
        line(q - Point((0, 0.4)), s + Point((0, -U - 1.0)))
    d.add(elm.Ground().at(s + Point((0, -U - 1.0))))
    return s

# ---- series branch ----
p1 = Point((0, 0))
d.add(elm.Dot(open=True).at(p1).label('p1', 'left'))
n1 = p1 + Point((1.5, 0))
line(p1, n1)
d.add(elm.Dot().at(n1))
a = n1 + Point((U, 0))
d.add(elm.Resistor().endpoints(n1, a).label('Rs'))
b = a + Point((U, 0))
d.add(elm.Inductor2(loops=3).endpoints(a, b).label('Ls'))
c = parallel_horizontal(b, elm.Inductor2(loops=2).label('Lskin1'), elm.Resistor().label('Rskin1', 'bottom'))
e = parallel_horizontal(c, elm.Inductor2(loops=2).label('Lskin2'), elm.Resistor().label('Rskin2', 'bottom'))
n2 = e + Point((0.5, 0))
line(e, n2)
d.add(elm.Dot().at(n2))
p2 = n2 + Point((1.5, 0))
line(n2, p2)
d.add(elm.Dot(open=True).at(p2).label('p2', 'right'))

# ---- Cs across the coil ----
h = 2.0
line(n1, n1 + Point((0, h)))
line(n2, n2 + Point((0, h)))
mid = (n1.x + n2.x) / 2
line(n1 + Point((0, h)), Point((mid - U / 2, h)))
d.add(elm.Capacitor().endpoints(Point((mid - U / 2, h)), Point((mid + U / 2, h))).label('Cs'))
line(Point((mid + U / 2, h)), n2 + Point((0, h)))

# ---- substrate networks ----
s1 = shunt_network(n1, 'Cox1', 'Rsi1', 'Csi1', 's1', 'left')
s2 = shunt_network(n2, 'Cox2', 'Rsi2', 'Csi2', 's2', 'right')

# ---- substrate coupling s1 - s2 (only used if needed) ----
w = U + 1.0
left = Point((mid - w / 2, s1.y))
line(s1 + Point((0.8, 0)), left)
right = parallel_horizontal(left, elm.Resistor().label('Rsub12'), elm.Capacitor().label('Csub12', 'bottom'), gap=1.4)
line(right, s2 - Point((0.8, 0)))

d.save(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'inductor_fit_model.png'), dpi=110, transparent=False)
