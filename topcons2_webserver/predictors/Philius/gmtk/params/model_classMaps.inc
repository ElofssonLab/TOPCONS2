%##############################################################################
%################## Copyright 2008 University of Washington ###################
%##############################################################################



    % -------------------------------------------------------------------------
    % each state is then mapped to two different 'classes' : one according
    % to its amino-acid emission probabilities (the 'emitClass'), and one
    % according to its duration distribution (the 'lengthClass')

    variable: emitClass {
	type: discrete hidden cardinality EMIT_CLASS_CARD;
	weight: scale WEIGHT_emitClass;
	switchingparents: nil;
	conditionalparents: state(0) using DeterministicCPT("state2emitClass");
    }

    variable: lengthClass {
	type: discrete hidden cardinality LENGTH_CLASS_CARD;
	weight: scale WEIGHT_lengthClass;
	switchingparents: nil;
	conditionalparents: state(0) using DeterministicCPT("state2lengthClass");
    }

    % -------------------------------------------------------------------------
    % the state is also mapped to a 'topoClass' which has slightly lower
    % cardinality than the 'stateLabel' and from this we also set the 'spBit'
    % binary variable ...

    variable: topoClass {
        type: discrete hidden cardinality TOPO_CLASS_CARD;
        switchingparents: nil;
        conditionalparents: state(0) using DeterministicCPT("state2topoClass");
    }

    % -------------------------------------------------------------------------
