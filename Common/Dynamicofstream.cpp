/*
 * Dynamicofstream.cpp
 *
 *  Created on: Jun 19, 2013
 *      Author: cjustin
 */

#include "Dynamicofstream.h"
#include "gzstream.h"
#include "StringUtil.h"
#include <iostream>
#include <fstream>

Dynamicofstream::Dynamicofstream(const string &filename)
{
	if (endsWith(filename, ".gz")) {
		filestream = new ogzstream(filename.c_str(), ios::out);
		gz = true;
	} else {
		filestream = new ofstream(filename.c_str(), ios::out);
		gz = false;
	}
	assert(filestream->good());
}

ostream& Dynamicofstream::operator <<(const string& o)
{
	*filestream << o;
	return *filestream;
}

ostream& Dynamicofstream::operator <<(unsigned o)
{
	*filestream << o;
	return *filestream;
}

void Dynamicofstream::close()
{
	filestream->flush();
	if (gz) {
		ogzstream *temp = dynamic_cast<ogzstream*>(filestream);
		temp->close();
	} else {
		ofstream *temp = dynamic_cast<ofstream*>(filestream);
		temp->close();
	}
}

Dynamicofstream::~Dynamicofstream()
{
	close();
	delete filestream;
}

