#define PY_SSIZE_T_CLEAN
#include <Python.h>

static PyModuleDef _shewchuk_module = {
    PyModuleDef_HEAD_INIT,
    .m_name = "shewchuk",
    .m_doc = PyDoc_STR("Robust floating point operations."),
    .m_size = -1,
};

PyMODINIT_FUNC PyInit__shewchuk(void) {
  return PyModule_Create(&_shewchuk_module);
}
