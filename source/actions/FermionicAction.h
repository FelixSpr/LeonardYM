#ifndef FERMIONICACTION_H_
#define FERMIONICACTION_H_

#include "Energy.h"
#include "hmc_forces/Force.h"
#include "dirac_operators/DiracOperator.h"

namespace Update {

class FermionicAction : public Update::Energy, public Update::Force {
public:
	FermionicAction(DiracOperator* _diracOperator);
	~FermionicAction();

	DiracOperator* getDiracOperator() const;
	void setDiracOperator(DiracOperator* _diracOperator);
protected:
	DiracOperator* diracOperator;
};

} /* namespace Update */
#endif /* FERMIONICACTION_H_ */
