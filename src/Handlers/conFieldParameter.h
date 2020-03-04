#ifndef CONFIELDPARAMETER_H
#define CONFIELDPARAMETER_H

#include "../CommonHandler.h"

#include "vHandler.h"
#include "Callback.h"
#include "Design.h"

class  conFieldParameter  : public Design  {
	int field_id;
	size_t Pars;
	int Par_size; ///< Parameter space dimension
	int *Par_sizes; ///< Parameter space dimensions on all the processors
	int *Par_disp; ///< Offsets in the Parameter vector for all the processors
	int mpi_size, mpi_rank;
	int CalculateNumberOfParameters ();
	bool FlagInDesignSpace(big_flag_t);
	bool InDesignSpace(size_t);
	int LocalParameters (int type, double * tab);
	big_flag_t flag_mask, flag_value;
public:
	static std::string xmlname;
	int Init ();
	int NumberOfParameters ();
	int Parameters (int type, double * tab);
};

#endif // CONFIELDPARAMETER_H
