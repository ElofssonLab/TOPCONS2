%##############################################################################
%################## Copyright 2008 University of Washington ###################
%##############################################################################

    % -------------------------------------------------------------------------
    % new in this frame is a constraint on the ending state: this is 
    % implemented using a binary 'observed child', always observed to have the 
    % value 1, and the constraint is enforced using a CPT which defines which
    % parents can result in this child having the value 1

    variable: stateConstraint {
	type: discrete observed value 1 cardinality BINARY_CARD;
	switchingparents: nil;
	conditionalparents: state(0) using DeterministicCPT("finalState");
    }

    % -------------------------------------------------------------------------

    variable: pType {
        type: discrete hidden cardinality PROTEIN_TYPE_CARD;
        switchingparents: nil;
        conditionalparents: spBit(0), iNumTM(0) using DeterministicCPT("setPtype");
    }

    % -------------------------------------------------------------------------
    % now we want to score the final topology class ... 
    % first, we will have a 'topoBit' that is linked to the last 'topoClass'
    % through an observed-child mechanism ... and meanwhile the 'topoBit'
    % is *scored* by having the spBit and the iNumTM as parents ...

    variable: topoBit {
	type: discrete hidden cardinality BINARY_CARD;
	weight: scale WEIGHT_topoBit;
	switchingparents: nil;
	conditionalparents: spBit(0), iNumTM(0) using DenseCPT("topoBit_spBitNumTM");
    }

    variable: oc2 {
	type: discrete observed value 1 cardinality BINARY_CARD;
	switchingparents: nil;
	conditionalparents: topoClass(0), topoBit(0) using DeterministicCPT("topoClassAndBit");
    }

    % -------------------------------------------------------------------------
