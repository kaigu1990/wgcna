## My way to use wgcna

* `expression.txt` is expression profiles from proteins or mRNAs
* `trait.txt` is trait matrix，but is also a optional file. If there is no trait data, you can use group data instead. First column is sample id, second column is sample group id.
	
		S-control-1	S-control
		S-control-2	S-control
		S-control-3	S-control
		S-NHPS-L-1	S-NHPS-L
		S-NHPS-L-2	S-NHPS-L
		S-NHPS-L-3	S-NHPS-L
		S-NHPS-H-1	S-NHPS-H
		S-NHPS-H-2	S-NHPS-H
		S-NHPS-H-3	S-NHPS-H
		F-control-1	F-control
		F-control-2	F-control
		F-control-3	F-control
		F-NHPS-L-1	F-NHPS-L
		F-NHPS-L-2	F-NHPS-L
		F-NHPS-L-3	F-NHPS-L
		F-NHPS-H-1	F-NHPS-H
		F-NHPS-H-2	F-NHPS-H
		F-NHPS-H-3	F-NHPS-H

* All result files and photos from wgcna analysis in Demo
* The default values of WGCNA networkType/TOMtype is unsigned
* If trait_form is equal to 1, you shoud have trait matrix, otherwise group data instead
