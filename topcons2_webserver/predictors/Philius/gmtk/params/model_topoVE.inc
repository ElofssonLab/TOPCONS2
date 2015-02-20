%##############################################################################
%################## Copyright 2008 University of Washington ###################
%##############################################################################


    % -------------------------------------------------------------------------
    % in this second-pass test model, we will use "virtual evidence" to 
    % encourage the topology to be the one that had the maximum posterior
    % probability according to the first pass using the junction-tree
    % algorithm ...

    variable: topoVE {
	type: discrete observed value 1 cardinality BINARY_CARD;
	weight: scale WEIGHT_topoVE;
	switchingparents: nil;
	conditionalparents: topo10(0) using VirtualEvidenceCPT("topo_VE");
    }

    % -------------------------------------------------------------------------
