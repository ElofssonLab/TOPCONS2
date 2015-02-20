%##############################################################################
%################## Copyright 2008 University of Washington ###################
%##############################################################################

    % -------------------------------------------------------------------------
    % aa : this variable is the observed amino acid (mapped to an integer)
    % and is chosen according to the current emitClass ...

    variable: aa {
	type: discrete observed 0:0 cardinality AA_CARD;
	weight: scale WEIGHT_aa;
	switchingparents: nil;
	conditionalparents: emitClass(0) using DenseCPT("aaEmitByClass_CPT");
    }

    % -------------------------------------------------------------------------
