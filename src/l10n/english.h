/*
 * english.h
 *
 *  Created on: 21 Dec 2013
 *      Author: jowr
 */

#ifndef ENGLISH_H_
#define ENGLISH_H_

namespace CoolProp {

const std::string ERR_NOT_A_TWO_PHASE_FLUID("This is not a two-phase fluid. Please select another fluid or avoid two-phase functions.");
const std::string ERR_NOT_A_TWO_PHASE_STATE("This is not a two-phase state, update state with a two-phase set of inputs");
const std::string ERR_NOT_COMPRESSIBLE("This function is invalid for incompressible fluids.");
const std::string ERR_NOT_A_TWO_PHASE_FUNCTION("This function is invalid in the two-phase region.");

} /* namespace CoolProp */
#endif /* ENGLISH_H_ */
