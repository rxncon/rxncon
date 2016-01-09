# distutils: sources = test/simulation/rule_based/cython_cpp/HelloWorld.cpp


cdef extern from "HelloWorld.h":
    cdef cppclass HelloWorld:
        HelloWorld() except +
        void say()


cdef class PyHelloWorld:
    cdef HelloWorld *thisptr

    def __cinit__(self):
        pass

    def __dealloc__(self):
        del self.thisptr

    def say(self):
        self.thisptr.say()
