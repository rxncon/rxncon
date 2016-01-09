# import io
# from util.stdout_redirector import stdout_redirector
#
# import test.simulation.rule_based.cython_cpp.hello_world as hel
#
#
# def test_hello_world_object_instantiates():
#     hello = hel.PyHelloWorld()
#
#     assert isinstance(hello, hel.PyHelloWorld)
#     assert hasattr(hello, 'say')
#
#
# def test_hello_world_output_correct():
#     hello = hel.PyHelloWorld()
#
#     buf = io.BytesIO()
#
#     with stdout_redirector(buf):
#         hello.say()
#
#     output = buf.getvalue()
#
#     print(output)
#
