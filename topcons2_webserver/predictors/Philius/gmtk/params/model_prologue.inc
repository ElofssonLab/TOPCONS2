%##############################################################################
%################## Copyright 2008 University of Washington ###################
%##############################################################################


    % -------------------------------------------------------------------------
    % changeState : in frame0, this variable is simply initialized to TRUE

    variable: changeState {
	type: discrete hidden cardinality BINARY_CARD;
	switchingparents: nil;
	conditionalparents: nil using DeterministicCPT("setBinaryTRUE");
    }

    % -------------------------------------------------------------------------
    % topoChange : in frame0, this variable is simply initialized to TRUE

    variable: topoChange {
	type: discrete hidden cardinality BINARY_CARD;
	switchingparents: nil;
	conditionalparents: nil using DeterministicCPT("setBinaryTRUE");
    }

    variable: topo10 {
	type: discrete hidden cardinality TOPO_10_CARD;
	switchingparents: nil;
	conditionalparents: 
		topoChange(0), topoClass(0) using DeterministicCPT("two2oneMap");
    }

    % -------------------------------------------------------------------------
    % state : this variable indicates what state we are currently in ... 
    % in frame0, the state will be chosen based on the probability of
    % _starting_ in a particular state

    variable: state {
	type: discrete hidden cardinality STATE_CARD;
	weight: scale WEIGHT_state;
	switchingparents: nil;
	conditionalparents: nil using DenseCPT("startState_CPT");
    }

    % -------------------------------------------------------------------------
    % the state is also mapped to a 'topoClass' which has slightly lower
    % cardinality than the 'stateLabel' and from this we also set the 'spBit'
    % and the 'iNumTM' variables from which we will later derive the 'pType'
    % in the epilogue ...

    variable: spBit {
        type: discrete hidden cardinality BINARY_CARD;
        switchingparents: nil;
        conditionalparents: topoClass(0) using DeterministicCPT("isSP");
    }

    variable: iNumTM {
        type: discrete hidden cardinality NUM_TM_CARD;
        switchingparents: nil;
        conditionalparents: nil using DeterministicCPT("setNumTMZero");
    }

    % -------------------------------------------------------------------------
    % stateCountDown : this variable is used to count *down* the duration
    % (aka length) of a particular state ... so the parents are :
    %	- from current  frame : changeState (as a switching parent)
    %		if changeState=FALSE, return a bad-value (since it should
    %			be impossible to get to this point)
    %		if changeState=TRUE,  choose duration based on lengthClass

    variable: stateCountDown {
	type: discrete hidden cardinality MAX_STATE_LENGTH_CARD;
	weight: scale WEIGHT_stateCountDown;
	switchingparents: changeState(0) using mapping ("oneD_map");
	conditionalparents: 
	    nil using DeterministicCPT("badLength")		% 0 
	  | lengthClass(0) using DenseCPT("length_CPT");	% 1 
    }

    % -------------------------------------------------------------------------
