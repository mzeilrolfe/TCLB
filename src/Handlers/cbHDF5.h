#ifndef CBHDF5_H
#define CBHDF5_H

#include "../CommonHandler.h"

#include "vHandler.h"
#include "Callback.h"

class  cbHDF5  : public  Callback  {
	std::string nm;
	name_set s;
	public:
	static std::string xmlname;
int Init ();
int DoIt ();
};

#endif // CBHDF5_H
