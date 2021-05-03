var searchIndex = JSON.parse('{\
"rustyread":{"doc":"","i":[[0,"alignment","rustyread","Function to align sequence",null,null],[5,"edit_distance","rustyread::alignment","",null,[[],["u64",15]]],[5,"align","","",null,[[]]],[0,"cli","rustyread","All stuff relate to command line",null,null],[0,"simulate","rustyread::cli","All stuff relate to simulate subcommand",null,null],[3,"Quantity","rustyread::cli::simulate","Store quantity as coverage of number of base",null,null],[11,"number_of_base","","Convert Quantity in a number of base, if base is set …",0,[[["u64",15]],["u64",15]]],[3,"Duo","","Store a pair of value, can be parse from str if it\'s …",null,null],[12,"0","","",1,null],[12,"1","","",1,null],[3,"Trio","","Store a trio of value, can be parse from str if it\'s …",null,null],[12,"0","","",2,null],[12,"1","","",2,null],[12,"2","","",2,null],[5,"found_model","","Found path to model file",null,[[["string",3]],[["pathbuf",3],["result",6]]]],[3,"Command","","Struct use to parse simulate subcommand argument",null,null],[12,"reference_path","","Path to reference sequence in fasta format",3,null],[12,"output_path","","Path to reference sequence in fasta format",3,null],[12,"quantity","","Quantity of base rustyread have to generate",3,null],[12,"length","","Read length distribution parameter",3,null],[12,"identity","","Identity distribution parameter",3,null],[12,"error_model","","Error model used",3,null],[12,"qscore_model","","Qualtity score model used",3,null],[12,"seed","","Seed used",3,null],[12,"start_adapter","","Start adapter parameter",3,null],[12,"end_adapter","","End adapter parameter",3,null],[12,"start_adapter_seq","","Start adapter sequence",3,null],[12,"end_adapter_seq","","End adapter sequence",3,null],[12,"junk","","Junk reads parameter",3,null],[12,"random","","Random reads parameter",3,null],[12,"chimera","","Chimera parameter",3,null],[12,"glitches","","Glitches parameter",3,null],[12,"small_plasmid_bias","","Small plasmid bias or not",3,null],[12,"nb_base_store","","Limit memory usage",3,null],[3,"Command","rustyread::cli","",null,null],[12,"subcmd","","Subcommand call",4,null],[12,"threads","","Number of thread rustyread can use",4,null],[12,"verbosity","","Verbosity level",4,null],[4,"SubCommand","","",null,null],[13,"Simulate","","",5,null],[5,"i82level","","Convert verbosity level (number of v) is log::Level",null,[[["i8",15]],[["level",4],["option",4]]]],[5,"set_nb_threads","","set number of global rayon thread pool",null,[[["usize",15]]]],[0,"error","rustyread","All stuff relate to error",null,null],[0,"cli","rustyread::error","Command line interface error",null,null],[4,"Cli","rustyread::error::cli","Enum to manage error polymorphism",null,null],[13,"CantParseQuantity","","quantity didn\'t match to pattern \\\\d+[KMGx]?",6,null],[13,"CantParseDuo","","Cant parse a duo of value",6,null],[13,"CantParseTrio","","Cant parse a trio of value",6,null],[13,"CantFoundModelPath","","Cant found model path",6,null],[13,"SmallPlasmidBias","","Small plasmid bias",6,null],[0,"model","rustyread::error","Model error",null,null],[4,"Model","rustyread::error::model","Enum to manage error polymorphism",null,null],[13,"ErrorParsing","","Error durring error model parsing",7,null],[13,"QualityParsing","","Error durring quality model parsing",7,null],[13,"QualityNotMinimalCigarString","","Quality model not contains minimal cigar string",7,null],[13,"QualityCigarLenNotOdd","","Cigar string length must be odd",7,null],[13,"LengthParamMustBeUpperThan0","","Length model parameter must be upper than 0.0",7,null],[13,"IdentityParamMustBeUpperThan0","","Identity model parameter must be upper than 0.0",7,null],[4,"Error","rustyread::error","Enum to manage error polymorphism",null,null],[13,"Cli","","Error related to command line interface",8,null],[13,"Model","","Error related to model",8,null],[0,"model","rustyread","Manage model",null,null],[0,"adapter","rustyread::model","Model to get sequence adapter",null,null],[3,"Adapter","rustyread::model::adapter","Struct to get sequence adapter",null,null],[11,"new","","Create model from parameter",9,[[["f64",15],["vec",3],["u8",15]],["result",6]]],[11,"get_start","","",9,[[],[["vec",3],["u8",15]]]],[11,"get_end","","",9,[[],[["vec",3],["u8",15]]]],[11,"max_len","","",9,[[],["usize",15]]],[0,"error","rustyread::model","Model to add error sequence",null,null],[3,"Error","rustyread::model::error","Struct to load and apply error model",null,null],[11,"from_stream","","Load model from an stdin",10,[[],["result",6]]],[11,"random","","Setup a random error model",10,[[["usize",15]]]],[11,"add_errors_to_kmer","","Add error to a kmer",10,[[]]],[11,"k","","Kmer length of model",10,[[],["usize",15]]],[5,"random_error","","Add a single random error in a kmer",null,[[],[["vec",3],["u8",15]]]],[0,"glitch","rustyread::model","Model to add glitches sequence",null,null],[3,"Glitch","rustyread::model::glitch","Struct to generate glitches in fragment",null,null],[11,"new","","Create model from parameter",11,[[["f64",15]],[["glitch",3],["result",6]]]],[11,"get_glitch","","Get glitch",11,[[],["option",4]]],[0,"identity","rustyread::model","Model to get read identity",null,null],[3,"Identity","rustyread::model::identity","Struct to generate length of fragment",null,null],[11,"new","","Create model from parameter",12,[[["f64",15]],[["identity",3],["result",6]]]],[11,"get_identity","","Get identity from model",12,[[],["f64",15]]],[0,"length","rustyread::model","Model to get length of reads",null,null],[3,"Length","rustyread::model::length","Struct to generate length of fragment",null,null],[11,"new","","Create model from parameter",13,[[["f64",15]],[["result",6],["length",3]]]],[11,"get_length","","Get length from model",13,[[],["u64",15]]],[0,"quality","rustyread::model","Model to generate qscore",null,null],[3,"Quality","rustyread::model::quality","Struct to load and apply quality model",null,null],[11,"from_stream","","Load model from an stdin",14,[[],["result",6]]],[11,"random","","Build a random quality score model",14,[[]]],[11,"ideal","","Build an ideal quality score model",14,[[]]],[11,"get_qscore","","Generate error associate to a cigar string with odd length",14,[[],[["result",6],["u8",15]]]],[11,"max_k","","Kmer length of model",14,[[],["usize",15]]],[0,"references","rustyread","A collections of sequence to store reference sequence",null,null],[3,"Reference","rustyread::references","Store a reference sequence",null,null],[12,"id","","",15,null],[12,"seq","","",15,null],[12,"revcomp","","",15,null],[12,"circular","","",15,null],[11,"new","","Build a new refrence",15,[[["box",3],["bool",15],["string",3]]]],[3,"References","","A collections of sequence",null,null],[12,"sequences","","",16,null],[12,"dist","","",16,null],[11,"from_stream","","Read a collection of sequence in fasta format from an …",16,[[],["result",6]]],[11,"from_stream_adjusted_weight","","Same as from_stream but small sequence have increase …",16,[[["bool",15],["length",3]],["result",6]]],[11,"choose_reference","","Randomly get a reference index and strand according to …",16,[[]]],[0,"simulate","rustyread","Simulate reads",null,null],[0,"description","rustyread::simulate","Manage read description",null,null],[4,"ReadType","rustyread::simulate::description","An enum to represent type of read",null,null],[13,"Real","","",17,null],[13,"Junk","","",17,null],[13,"Random","","",17,null],[3,"Origin","","Store information about origin of read",null,null],[12,"ref_id","","",18,null],[12,"strand","","",18,null],[12,"start","","",18,null],[12,"end","","",18,null],[12,"read_type","","",18,null],[11,"reference","","",18,[[["usize",15],["string",3],["char",15]]]],[11,"junk","","",18,[[["usize",15]]]],[11,"random","","",18,[[["usize",15]]]],[3,"Description","","Store information about read",null,null],[12,"origin","","",19,null],[12,"chimera","","",19,null],[12,"length","","",19,null],[12,"identity","","",19,null],[11,"new","","",19,[[["usize",15],["f64",15],["origin",3],["option",4]]]],[0,"error","rustyread::simulate","Add error on reads",null,null],[6,"Seq","rustyread::simulate::error","",null,null],[6,"Cigar","","",null,null],[3,"Change","","A struct to represent a change in a sequence",null,null],[11,"new","","Build a new Change warning CIGAR String and edit_distance …",20,[[["usize",15],["seq",6]]]],[11,"from_seq","","Build a new Change from a begin, end, a change sequence …",20,[[["usize",15],["seq",6]]]],[11,"begin","","Get begin of Change",20,[[],["usize",15]]],[11,"end_raw","","Get end of change in raw sequence",20,[[],["usize",15]]],[11,"end_err","","Get end of change in erroneous sequence",20,[[],["usize",15]]],[11,"seq","","Get erroneous sequence",20,[[],["seq",6]]],[11,"cigar","","Get cigar string",20,[[],["cigar",6]]],[11,"edit","","Get edit distance",20,[[],["f64",15]]],[11,"contain","","Return true if other is contain in self This function …",20,[[["change",3]],["bool",15]]],[11,"overlap","","Return True if self and Change overlap This function …",20,[[["change",3]],["bool",15]]],[11,"merge","","Try to merge two Change. This function assume other.begin …",20,[[["change",3]],["f64",15]]],[11,"update_align","","Update alignment info of error",20,[[]]],[6,"Changes","","",null,null],[5,"number_of_edit","","From identity an seq length compute number of edit",null,[[["f64",15],["usize",15]],["f64",15]]],[5,"sequence","","Apply error on read",null,[[["f64",15],["error",3],["glitch",3]]]],[5,"add_glitches","","Create Change correspond to glitches",null,[[["glitch",3],["changes",6]]]],[5,"add_error","","Add Change correspond to error in changes",null,[[["usize",15],["f64",15],["error",3],["changes",6]]]],[0,"fragments","rustyread::simulate","Generate fragments",null,null],[3,"Fragments","rustyread::simulate::fragments","An iterator produce fragment, a description and a seed",null,null],[11,"new","","Create a new Fragments",21,[[["length",3],["identity",3],["u64",15],["references",3]]]],[11,"get_read_type","","Get the read type",21,[[],["readtype",4]]],[11,"is_chimera","","Return true fragment must be a chimera",21,[[],["bool",15]]],[11,"generate_fragment","","Produce a fragment",21,[[]]],[0,"quality","rustyread::simulate","Generate quality for simulate",null,null],[5,"generate_quality","rustyread::simulate::quality","Generate quality string",null,[[["quality",3]],[["result",6],["vec",3]]]],[5,"simulate","rustyread::simulate","main simulate function",null,[[["command",3]],["result",6]]],[5,"random_base","rustyread","Get a random base",null,[[],["u8",15]]],[5,"random_base_diff","","Get a random base diffrent than nuc",null,[[["u8",15]],["u8",15]]],[5,"random_seq","","Get random sequences",null,[[["usize",15]],[["vec",3],["u8",15]]]],[11,"from","rustyread::cli::simulate","",0,[[]]],[11,"into","","",0,[[]]],[11,"borrow","","",0,[[]]],[11,"borrow_mut","","",0,[[]]],[11,"try_from","","",0,[[],["result",4]]],[11,"try_into","","",0,[[],["result",4]]],[11,"type_id","","",0,[[],["typeid",3]]],[11,"vzip","","",0,[[]]],[11,"init","","",0,[[],["usize",15]]],[11,"deref","","",0,[[["usize",15]]]],[11,"deref_mut","","",0,[[["usize",15]]]],[11,"drop","","",0,[[["usize",15]]]],[11,"to_subset","","",0,[[],["option",4]]],[11,"is_in_subset","","",0,[[],["bool",15]]],[11,"to_subset_unchecked","","",0,[[]]],[11,"from_subset","","",0,[[]]],[11,"from","","",1,[[]]],[11,"into","","",1,[[]]],[11,"borrow","","",1,[[]]],[11,"borrow_mut","","",1,[[]]],[11,"try_from","","",1,[[],["result",4]]],[11,"try_into","","",1,[[],["result",4]]],[11,"type_id","","",1,[[],["typeid",3]]],[11,"vzip","","",1,[[]]],[11,"init","","",1,[[],["usize",15]]],[11,"deref","","",1,[[["usize",15]]]],[11,"deref_mut","","",1,[[["usize",15]]]],[11,"drop","","",1,[[["usize",15]]]],[11,"to_subset","","",1,[[],["option",4]]],[11,"is_in_subset","","",1,[[],["bool",15]]],[11,"to_subset_unchecked","","",1,[[]]],[11,"from_subset","","",1,[[]]],[11,"from","","",2,[[]]],[11,"into","","",2,[[]]],[11,"borrow","","",2,[[]]],[11,"borrow_mut","","",2,[[]]],[11,"try_from","","",2,[[],["result",4]]],[11,"try_into","","",2,[[],["result",4]]],[11,"type_id","","",2,[[],["typeid",3]]],[11,"vzip","","",2,[[]]],[11,"init","","",2,[[],["usize",15]]],[11,"deref","","",2,[[["usize",15]]]],[11,"deref_mut","","",2,[[["usize",15]]]],[11,"drop","","",2,[[["usize",15]]]],[11,"to_subset","","",2,[[],["option",4]]],[11,"is_in_subset","","",2,[[],["bool",15]]],[11,"to_subset_unchecked","","",2,[[]]],[11,"from_subset","","",2,[[]]],[11,"from","","",3,[[]]],[11,"into","","",3,[[]]],[11,"borrow","","",3,[[]]],[11,"borrow_mut","","",3,[[]]],[11,"try_from","","",3,[[],["result",4]]],[11,"try_into","","",3,[[],["result",4]]],[11,"type_id","","",3,[[],["typeid",3]]],[11,"vzip","","",3,[[]]],[11,"init","","",3,[[],["usize",15]]],[11,"deref","","",3,[[["usize",15]]]],[11,"deref_mut","","",3,[[["usize",15]]]],[11,"drop","","",3,[[["usize",15]]]],[11,"to_subset","","",3,[[],["option",4]]],[11,"is_in_subset","","",3,[[],["bool",15]]],[11,"to_subset_unchecked","","",3,[[]]],[11,"from_subset","","",3,[[]]],[11,"from","rustyread::cli","",4,[[]]],[11,"into","","",4,[[]]],[11,"borrow","","",4,[[]]],[11,"borrow_mut","","",4,[[]]],[11,"try_from","","",4,[[],["result",4]]],[11,"try_into","","",4,[[],["result",4]]],[11,"type_id","","",4,[[],["typeid",3]]],[11,"vzip","","",4,[[]]],[11,"init","","",4,[[],["usize",15]]],[11,"deref","","",4,[[["usize",15]]]],[11,"deref_mut","","",4,[[["usize",15]]]],[11,"drop","","",4,[[["usize",15]]]],[11,"to_subset","","",4,[[],["option",4]]],[11,"is_in_subset","","",4,[[],["bool",15]]],[11,"to_subset_unchecked","","",4,[[]]],[11,"from_subset","","",4,[[]]],[11,"from","","",5,[[]]],[11,"into","","",5,[[]]],[11,"borrow","","",5,[[]]],[11,"borrow_mut","","",5,[[]]],[11,"try_from","","",5,[[],["result",4]]],[11,"try_into","","",5,[[],["result",4]]],[11,"type_id","","",5,[[],["typeid",3]]],[11,"vzip","","",5,[[]]],[11,"init","","",5,[[],["usize",15]]],[11,"deref","","",5,[[["usize",15]]]],[11,"deref_mut","","",5,[[["usize",15]]]],[11,"drop","","",5,[[["usize",15]]]],[11,"to_subset","","",5,[[],["option",4]]],[11,"is_in_subset","","",5,[[],["bool",15]]],[11,"to_subset_unchecked","","",5,[[]]],[11,"from_subset","","",5,[[]]],[11,"from","rustyread::error::cli","",6,[[]]],[11,"into","","",6,[[]]],[11,"to_string","","",6,[[],["string",3]]],[11,"borrow","","",6,[[]]],[11,"borrow_mut","","",6,[[]]],[11,"try_from","","",6,[[],["result",4]]],[11,"try_into","","",6,[[],["result",4]]],[11,"type_id","","",6,[[],["typeid",3]]],[11,"vzip","","",6,[[]]],[11,"init","","",6,[[],["usize",15]]],[11,"deref","","",6,[[["usize",15]]]],[11,"deref_mut","","",6,[[["usize",15]]]],[11,"drop","","",6,[[["usize",15]]]],[11,"to_subset","","",6,[[],["option",4]]],[11,"is_in_subset","","",6,[[],["bool",15]]],[11,"to_subset_unchecked","","",6,[[]]],[11,"from_subset","","",6,[[]]],[11,"from","rustyread::error::model","",7,[[]]],[11,"into","","",7,[[]]],[11,"to_string","","",7,[[],["string",3]]],[11,"borrow","","",7,[[]]],[11,"borrow_mut","","",7,[[]]],[11,"try_from","","",7,[[],["result",4]]],[11,"try_into","","",7,[[],["result",4]]],[11,"type_id","","",7,[[],["typeid",3]]],[11,"vzip","","",7,[[]]],[11,"init","","",7,[[],["usize",15]]],[11,"deref","","",7,[[["usize",15]]]],[11,"deref_mut","","",7,[[["usize",15]]]],[11,"drop","","",7,[[["usize",15]]]],[11,"to_subset","","",7,[[],["option",4]]],[11,"is_in_subset","","",7,[[],["bool",15]]],[11,"to_subset_unchecked","","",7,[[]]],[11,"from_subset","","",7,[[]]],[11,"from","rustyread::error","",8,[[]]],[11,"into","","",8,[[]]],[11,"to_string","","",8,[[],["string",3]]],[11,"borrow","","",8,[[]]],[11,"borrow_mut","","",8,[[]]],[11,"try_from","","",8,[[],["result",4]]],[11,"try_into","","",8,[[],["result",4]]],[11,"type_id","","",8,[[],["typeid",3]]],[11,"vzip","","",8,[[]]],[11,"init","","",8,[[],["usize",15]]],[11,"deref","","",8,[[["usize",15]]]],[11,"deref_mut","","",8,[[["usize",15]]]],[11,"drop","","",8,[[["usize",15]]]],[11,"to_subset","","",8,[[],["option",4]]],[11,"is_in_subset","","",8,[[],["bool",15]]],[11,"to_subset_unchecked","","",8,[[]]],[11,"from_subset","","",8,[[]]],[11,"from","rustyread::model::adapter","",9,[[]]],[11,"into","","",9,[[]]],[11,"borrow","","",9,[[]]],[11,"borrow_mut","","",9,[[]]],[11,"try_from","","",9,[[],["result",4]]],[11,"try_into","","",9,[[],["result",4]]],[11,"type_id","","",9,[[],["typeid",3]]],[11,"vzip","","",9,[[]]],[11,"init","","",9,[[],["usize",15]]],[11,"deref","","",9,[[["usize",15]]]],[11,"deref_mut","","",9,[[["usize",15]]]],[11,"drop","","",9,[[["usize",15]]]],[11,"to_subset","","",9,[[],["option",4]]],[11,"is_in_subset","","",9,[[],["bool",15]]],[11,"to_subset_unchecked","","",9,[[]]],[11,"from_subset","","",9,[[]]],[11,"from","rustyread::model::error","",10,[[]]],[11,"into","","",10,[[]]],[11,"borrow","","",10,[[]]],[11,"borrow_mut","","",10,[[]]],[11,"try_from","","",10,[[],["result",4]]],[11,"try_into","","",10,[[],["result",4]]],[11,"type_id","","",10,[[],["typeid",3]]],[11,"vzip","","",10,[[]]],[11,"init","","",10,[[],["usize",15]]],[11,"deref","","",10,[[["usize",15]]]],[11,"deref_mut","","",10,[[["usize",15]]]],[11,"drop","","",10,[[["usize",15]]]],[11,"to_subset","","",10,[[],["option",4]]],[11,"is_in_subset","","",10,[[],["bool",15]]],[11,"to_subset_unchecked","","",10,[[]]],[11,"from_subset","","",10,[[]]],[11,"from","rustyread::model::glitch","",11,[[]]],[11,"into","","",11,[[]]],[11,"borrow","","",11,[[]]],[11,"borrow_mut","","",11,[[]]],[11,"try_from","","",11,[[],["result",4]]],[11,"try_into","","",11,[[],["result",4]]],[11,"type_id","","",11,[[],["typeid",3]]],[11,"vzip","","",11,[[]]],[11,"init","","",11,[[],["usize",15]]],[11,"deref","","",11,[[["usize",15]]]],[11,"deref_mut","","",11,[[["usize",15]]]],[11,"drop","","",11,[[["usize",15]]]],[11,"to_subset","","",11,[[],["option",4]]],[11,"is_in_subset","","",11,[[],["bool",15]]],[11,"to_subset_unchecked","","",11,[[]]],[11,"from_subset","","",11,[[]]],[11,"from","rustyread::model::identity","",12,[[]]],[11,"into","","",12,[[]]],[11,"borrow","","",12,[[]]],[11,"borrow_mut","","",12,[[]]],[11,"try_from","","",12,[[],["result",4]]],[11,"try_into","","",12,[[],["result",4]]],[11,"type_id","","",12,[[],["typeid",3]]],[11,"vzip","","",12,[[]]],[11,"init","","",12,[[],["usize",15]]],[11,"deref","","",12,[[["usize",15]]]],[11,"deref_mut","","",12,[[["usize",15]]]],[11,"drop","","",12,[[["usize",15]]]],[11,"to_subset","","",12,[[],["option",4]]],[11,"is_in_subset","","",12,[[],["bool",15]]],[11,"to_subset_unchecked","","",12,[[]]],[11,"from_subset","","",12,[[]]],[11,"from","rustyread::model::length","",13,[[]]],[11,"into","","",13,[[]]],[11,"borrow","","",13,[[]]],[11,"borrow_mut","","",13,[[]]],[11,"try_from","","",13,[[],["result",4]]],[11,"try_into","","",13,[[],["result",4]]],[11,"type_id","","",13,[[],["typeid",3]]],[11,"vzip","","",13,[[]]],[11,"init","","",13,[[],["usize",15]]],[11,"deref","","",13,[[["usize",15]]]],[11,"deref_mut","","",13,[[["usize",15]]]],[11,"drop","","",13,[[["usize",15]]]],[11,"to_subset","","",13,[[],["option",4]]],[11,"is_in_subset","","",13,[[],["bool",15]]],[11,"to_subset_unchecked","","",13,[[]]],[11,"from_subset","","",13,[[]]],[11,"from","rustyread::model::quality","",14,[[]]],[11,"into","","",14,[[]]],[11,"borrow","","",14,[[]]],[11,"borrow_mut","","",14,[[]]],[11,"try_from","","",14,[[],["result",4]]],[11,"try_into","","",14,[[],["result",4]]],[11,"type_id","","",14,[[],["typeid",3]]],[11,"vzip","","",14,[[]]],[11,"init","","",14,[[],["usize",15]]],[11,"deref","","",14,[[["usize",15]]]],[11,"deref_mut","","",14,[[["usize",15]]]],[11,"drop","","",14,[[["usize",15]]]],[11,"to_subset","","",14,[[],["option",4]]],[11,"is_in_subset","","",14,[[],["bool",15]]],[11,"to_subset_unchecked","","",14,[[]]],[11,"from_subset","","",14,[[]]],[11,"from","rustyread::references","",15,[[]]],[11,"into","","",15,[[]]],[11,"borrow","","",15,[[]]],[11,"borrow_mut","","",15,[[]]],[11,"try_from","","",15,[[],["result",4]]],[11,"try_into","","",15,[[],["result",4]]],[11,"type_id","","",15,[[],["typeid",3]]],[11,"vzip","","",15,[[]]],[11,"init","","",15,[[],["usize",15]]],[11,"deref","","",15,[[["usize",15]]]],[11,"deref_mut","","",15,[[["usize",15]]]],[11,"drop","","",15,[[["usize",15]]]],[11,"to_subset","","",15,[[],["option",4]]],[11,"is_in_subset","","",15,[[],["bool",15]]],[11,"to_subset_unchecked","","",15,[[]]],[11,"from_subset","","",15,[[]]],[11,"from","","",16,[[]]],[11,"into","","",16,[[]]],[11,"borrow","","",16,[[]]],[11,"borrow_mut","","",16,[[]]],[11,"try_from","","",16,[[],["result",4]]],[11,"try_into","","",16,[[],["result",4]]],[11,"type_id","","",16,[[],["typeid",3]]],[11,"vzip","","",16,[[]]],[11,"init","","",16,[[],["usize",15]]],[11,"deref","","",16,[[["usize",15]]]],[11,"deref_mut","","",16,[[["usize",15]]]],[11,"drop","","",16,[[["usize",15]]]],[11,"to_subset","","",16,[[],["option",4]]],[11,"is_in_subset","","",16,[[],["bool",15]]],[11,"to_subset_unchecked","","",16,[[]]],[11,"from_subset","","",16,[[]]],[11,"from","rustyread::simulate::description","",17,[[]]],[11,"into","","",17,[[]]],[11,"to_owned","","",17,[[]]],[11,"clone_into","","",17,[[]]],[11,"borrow","","",17,[[]]],[11,"borrow_mut","","",17,[[]]],[11,"try_from","","",17,[[],["result",4]]],[11,"try_into","","",17,[[],["result",4]]],[11,"type_id","","",17,[[],["typeid",3]]],[11,"vzip","","",17,[[]]],[11,"init","","",17,[[],["usize",15]]],[11,"deref","","",17,[[["usize",15]]]],[11,"deref_mut","","",17,[[["usize",15]]]],[11,"drop","","",17,[[["usize",15]]]],[11,"to_subset","","",17,[[],["option",4]]],[11,"is_in_subset","","",17,[[],["bool",15]]],[11,"to_subset_unchecked","","",17,[[]]],[11,"from_subset","","",17,[[]]],[11,"from","","",18,[[]]],[11,"into","","",18,[[]]],[11,"to_owned","","",18,[[]]],[11,"clone_into","","",18,[[]]],[11,"to_string","","",18,[[],["string",3]]],[11,"borrow","","",18,[[]]],[11,"borrow_mut","","",18,[[]]],[11,"try_from","","",18,[[],["result",4]]],[11,"try_into","","",18,[[],["result",4]]],[11,"type_id","","",18,[[],["typeid",3]]],[11,"vzip","","",18,[[]]],[11,"init","","",18,[[],["usize",15]]],[11,"deref","","",18,[[["usize",15]]]],[11,"deref_mut","","",18,[[["usize",15]]]],[11,"drop","","",18,[[["usize",15]]]],[11,"to_subset","","",18,[[],["option",4]]],[11,"is_in_subset","","",18,[[],["bool",15]]],[11,"to_subset_unchecked","","",18,[[]]],[11,"from_subset","","",18,[[]]],[11,"from","","",19,[[]]],[11,"into","","",19,[[]]],[11,"to_owned","","",19,[[]]],[11,"clone_into","","",19,[[]]],[11,"to_string","","",19,[[],["string",3]]],[11,"borrow","","",19,[[]]],[11,"borrow_mut","","",19,[[]]],[11,"try_from","","",19,[[],["result",4]]],[11,"try_into","","",19,[[],["result",4]]],[11,"type_id","","",19,[[],["typeid",3]]],[11,"vzip","","",19,[[]]],[11,"init","","",19,[[],["usize",15]]],[11,"deref","","",19,[[["usize",15]]]],[11,"deref_mut","","",19,[[["usize",15]]]],[11,"drop","","",19,[[["usize",15]]]],[11,"to_subset","","",19,[[],["option",4]]],[11,"is_in_subset","","",19,[[],["bool",15]]],[11,"to_subset_unchecked","","",19,[[]]],[11,"from_subset","","",19,[[]]],[11,"from","rustyread::simulate::error","",20,[[]]],[11,"into","","",20,[[]]],[11,"to_owned","","",20,[[]]],[11,"clone_into","","",20,[[]]],[11,"borrow","","",20,[[]]],[11,"borrow_mut","","",20,[[]]],[11,"try_from","","",20,[[],["result",4]]],[11,"try_into","","",20,[[],["result",4]]],[11,"type_id","","",20,[[],["typeid",3]]],[11,"vzip","","",20,[[]]],[11,"init","","",20,[[],["usize",15]]],[11,"deref","","",20,[[["usize",15]]]],[11,"deref_mut","","",20,[[["usize",15]]]],[11,"drop","","",20,[[["usize",15]]]],[11,"to_subset","","",20,[[],["option",4]]],[11,"is_in_subset","","",20,[[],["bool",15]]],[11,"to_subset_unchecked","","",20,[[]]],[11,"from_subset","","",20,[[]]],[11,"from","rustyread::simulate::fragments","",21,[[]]],[11,"into","","",21,[[]]],[11,"into_iter","","",21,[[]]],[11,"borrow","","",21,[[]]],[11,"borrow_mut","","",21,[[]]],[11,"try_from","","",21,[[],["result",4]]],[11,"try_into","","",21,[[],["result",4]]],[11,"type_id","","",21,[[],["typeid",3]]],[11,"vzip","","",21,[[]]],[11,"init","","",21,[[],["usize",15]]],[11,"deref","","",21,[[["usize",15]]]],[11,"deref_mut","","",21,[[["usize",15]]]],[11,"drop","","",21,[[["usize",15]]]],[11,"to_subset","","",21,[[],["option",4]]],[11,"is_in_subset","","",21,[[],["bool",15]]],[11,"to_subset_unchecked","","",21,[[]]],[11,"from_subset","","",21,[[]]],[11,"from","rustyread::error","",8,[[["cli",4]]]],[11,"from","","",8,[[["model",4]]]],[11,"next","rustyread::simulate::fragments","",21,[[],["option",4]]],[11,"clone","rustyread::simulate::description","",17,[[],["readtype",4]]],[11,"clone","","",18,[[],["origin",3]]],[11,"clone","","",19,[[],["description",3]]],[11,"clone","rustyread::simulate::error","",20,[[],["change",3]]],[11,"eq","rustyread::cli::simulate","",0,[[["quantity",3]],["bool",15]]],[11,"ne","","",0,[[["quantity",3]],["bool",15]]],[11,"eq","","",1,[[["duo",3]],["bool",15]]],[11,"ne","","",1,[[["duo",3]],["bool",15]]],[11,"eq","","",2,[[["trio",3]],["bool",15]]],[11,"ne","","",2,[[["trio",3]],["bool",15]]],[11,"eq","rustyread::references","",15,[[["reference",3]],["bool",15]]],[11,"ne","","",15,[[["reference",3]],["bool",15]]],[11,"eq","rustyread::simulate::description","",17,[[["readtype",4]],["bool",15]]],[11,"eq","","",18,[[["origin",3]],["bool",15]]],[11,"ne","","",18,[[["origin",3]],["bool",15]]],[11,"eq","","",19,[[["description",3]],["bool",15]]],[11,"ne","","",19,[[["description",3]],["bool",15]]],[11,"eq","rustyread::simulate::error","",20,[[["change",3]],["bool",15]]],[11,"ne","","",20,[[["change",3]],["bool",15]]],[11,"fmt","rustyread::cli::simulate","",0,[[["formatter",3]],["result",6]]],[11,"fmt","","",1,[[["formatter",3]],["result",6]]],[11,"fmt","","",2,[[["formatter",3]],["result",6]]],[11,"fmt","","",3,[[["formatter",3]],["result",6]]],[11,"fmt","rustyread::cli","",4,[[["formatter",3]],["result",6]]],[11,"fmt","","",5,[[["formatter",3]],["result",6]]],[11,"fmt","rustyread::error::cli","",6,[[["formatter",3]],["result",6]]],[11,"fmt","rustyread::error::model","",7,[[["formatter",3]],["result",6]]],[11,"fmt","rustyread::error","",8,[[["formatter",3]],["result",6]]],[11,"fmt","rustyread::model::identity","",12,[[["formatter",3]],["result",6]]],[11,"fmt","rustyread::references","",15,[[["formatter",3]],["result",6]]],[11,"fmt","rustyread::simulate::description","",17,[[["formatter",3]],["result",6]]],[11,"fmt","","",18,[[["formatter",3]],["result",6]]],[11,"fmt","","",19,[[["formatter",3]],["result",6]]],[11,"fmt","rustyread::simulate::error","",20,[[["formatter",3]],["result",6]]],[11,"fmt","rustyread::error::cli","",6,[[["formatter",3]],["result",6]]],[11,"fmt","rustyread::error::model","",7,[[["formatter",3]],["result",6]]],[11,"fmt","rustyread::error","",8,[[["formatter",3]],["result",6]]],[11,"fmt","rustyread::simulate::description","",18,[[["formatter",3]],["result",6]]],[11,"fmt","","",19,[[["formatter",3]],["result",6]]],[11,"from_str","rustyread::cli::simulate","",0,[[["str",15]],["result",4]]],[11,"from_str","","",1,[[["str",15]],["result",4]]],[11,"from_str","","",2,[[["str",15]],["result",4]]],[11,"source","rustyread::error","",8,[[],[["option",4],["error",8]]]],[11,"into_app","rustyread::cli::simulate","",3,[[],["app",3]]],[11,"augment_clap","","",3,[[["app",3]],["app",3]]],[11,"into_app","rustyread::cli","",4,[[],["app",3]]],[11,"augment_clap","","",4,[[["app",3]],["app",3]]],[11,"into_app","","",5,[[],["app",3]]],[11,"augment_clap","","",5,[[["app",3]],["app",3]]],[11,"from_arg_matches","rustyread::cli::simulate","",3,[[["argmatches",3]]]],[11,"from_arg_matches","rustyread::cli","",4,[[["argmatches",3]]]],[11,"from_arg_matches","","",5,[[["argmatches",3]]]],[11,"augment_subcommands","","",5,[[["app",3]],["app",3]]],[11,"from_subcommand","","",5,[[["option",4]],["option",4]]]],"p":[[3,"Quantity"],[3,"Duo"],[3,"Trio"],[3,"Command"],[3,"Command"],[4,"SubCommand"],[4,"Cli"],[4,"Model"],[4,"Error"],[3,"Adapter"],[3,"Error"],[3,"Glitch"],[3,"Identity"],[3,"Length"],[3,"Quality"],[3,"Reference"],[3,"References"],[4,"ReadType"],[3,"Origin"],[3,"Description"],[3,"Change"],[3,"Fragments"]]}\
}');
addSearchOptions(searchIndex);initSearch(searchIndex);