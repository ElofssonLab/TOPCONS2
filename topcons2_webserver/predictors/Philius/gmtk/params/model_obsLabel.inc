%##############################################################################
%################## Copyright 2008 University of Washington ###################
%##############################################################################

    % -------------------------------------------------------------------------
    % stateLabel : the label is observed in the 2nd column of the input file

    variable: stateLabel {
	type: discrete observed 1:1 cardinality LABEL_CARD;
	switchingparents: nil;
	conditionalparents: nil using DenseCPT("internal:UnityScore");
    }

    % now we use an 'observed child' in order to implement the relationship
    % between 'state' and 'stateLabel'

    variable: oc {
	type: discrete observed value 1 cardinality BINARY_CARD;
	switchingparents: nil;
	conditionalparents: state(0), stateLabel(0) using DeterministicCPT("stateANDlabel");
    }

    % -------------------------------------------------------------------------
