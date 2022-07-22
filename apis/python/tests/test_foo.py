import tiledbsc.v1.general_utilities as gu

# TEMP TEMP TYPEGUARD EXPERIMENT
def test_foo():
    a = gu.A()
    b = gu.B()
    fa = gu.f(a)
    fb = gu.f(b)

    assert fa == fb
