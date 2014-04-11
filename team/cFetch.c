#include <stdio.h>
#include <stdlib.h>
#include <Python.h>
#include "methStruct.h"

//http://www.codingunit.com/c-tutorial-binary-file-io

void fetch(char* fileName, int seekStart, int blocks, 
	int minCov, int maxCov, PyObject* outMeth, PyObject* outContext) {
	FILE* f=fopen(fileName,"rb");
	int fseekError = fseek(f, (long int) seekStart, SEEK_SET);
	// Check for errors in fseek
	if(fseekError) {
		printf("Seek error\n");
	}
	methStruct* structArray = malloc(sizeof(methStruct)*blocks);
	size_t readRet = fread(structArray, sizeof(methStruct), 
		(size_t) blocks, f);
	if(readRet != blocks) {
		printf("Read failed.\n");
	}
	fclose(f);
	int i;
	double tmp;
	for(i=0; i<blocks; i++) {
		// Check for no coverage
		if(structArray[i].c == (unsigned short) 65535 && structArray[i].ct == (unsigned short) 65535) {
			PyList_SetItem(outMeth, (Py_ssize_t) i, Py_BuildValue("d", -1.0));
		// Check for not enough coverage
		} else if(structArray[i].ct < minCov) {
			PyList_SetItem(outMeth, (Py_ssize_t) i, Py_BuildValue("d", -1.0));
		// Check for too much coverage (from repeats)
		} else if(structArray[i].ct > maxCov) {
			PyList_SetItem(outMeth, (Py_ssize_t) i, Py_BuildValue("d", -1.0));
		} else {
			tmp = (double) structArray[i].c / (double) structArray[i].ct;
			PyList_SetItem(outMeth, (Py_ssize_t) i, Py_BuildValue("d", tmp));
		}
		PyList_SetItem(outContext, (Py_ssize_t) i, Py_BuildValue("i",(int) structArray[i].con));
	}
	free(structArray);
}

static PyObject* cFetch(PyObject* self, PyObject* args) {
	// Pass back a tuple of arrays
	char* fileName = NULL;
	int seekStart = 0;
	int blocks = 0;
	int minCov = 0;
	int maxCov = 0;
	if (!PyArg_ParseTuple(args, "siiii", &fileName, &seekStart, &blocks,
		&minCov, &maxCov)) {
		return NULL;
	}
	PyObject* outMeth = PyList_New((Py_ssize_t) blocks);
	PyObject* outContext = PyList_New((Py_ssize_t) blocks);
	fetch(fileName, seekStart, blocks, minCov, maxCov, outMeth, outContext);
	return Py_BuildValue("(OO)", outMeth, outContext);
}

static PyMethodDef FetchMethods[] = {
	{"cFetch", cFetch, METH_VARARGS, "Implements fetch in C"},
	{NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC
initcFetch(void) {
	(void) Py_InitModule("cFetch", FetchMethods);
}
