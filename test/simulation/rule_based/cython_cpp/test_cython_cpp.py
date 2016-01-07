import test.simulation.rule_based.cython_cpp.hello_world as hel


def test_hello_world_object_instantiates():
    hello = hel.PyHelloWorld()

    assert isinstance(hello, hel.PyHelloWorld)
    assert hasattr(hello, 'say')
