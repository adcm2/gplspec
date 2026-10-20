# Stage 1C candidate content manifest

- Base / Stage 1C parent: `b68bb219e94eac50ecc145c9eb6388ac5fb61f57`.
- Candidate at generation: working-tree changes from `b68bb219e94eac50ecc145c9eb6388ac5fb61f57`; no Stage 1C commit exists.
- Scope: all tracked GPLSpec production files, example sources/build inputs, test sources/scripts/build inputs, CMake modules, all tracked CI workflows, and the extracted tutorial source used by the build. Documentation, this manifest, `implementation_status.md`, and the consolidated review record are excluded to avoid recursive/self or review-only entries.
- Each row gives SHA-256 at the base commit and in the candidate. `UNTRACKED` means the candidate file is new and has no base blob.

| Base SHA-256 | Candidate SHA-256 | Path |
| --- | --- | --- |
| `UNTRACKED` | `1d9f616e2f787e4ff04797132a28063c38aa2aff10af3aa420b46aaa92bd0f4e` | `.github/workflows/cxx.yml` |
| `52d9d998819b96bfb546d4496ac592ac7ba7acc6abdc2f11c643f5c06ba64da9` | `52d9d998819b96bfb546d4496ac592ac7ba7acc6abdc2f11c643f5c06ba64da9` | `.github/workflows/pages.yml` |
| `71f294d6215c630697f01c919dac5f2fc7f00260491c741eb414ae6cf2cf78f1` | `71f294d6215c630697f01c919dac5f2fc7f00260491c741eb414ae6cf2cf78f1` | `CMakeLists.txt` |
| `a01a8aee615b1d63c2d536865b745c6411efeb1a666bc31119798c67c88d23ab` | `a01a8aee615b1d63c2d536865b745c6411efeb1a666bc31119798c67c88d23ab` | `_tutorials/tutorial5_phobos.md` |
| `677d2576b297e6e90ae633346403b05576c5dec98f33c83507deb555e06e17ce` | `677d2576b297e6e90ae633346403b05576c5dec98f33c83507deb555e06e17ce` | `cmake/ApplyPinnedOdrPatch.cmake` |
| `639237850da14e78fc0f7b59c62e05962c8373af67c1c678b2d8a95e0ddafb35` | `639237850da14e78fc0f7b59c62e05962c8373af67c1c678b2d8a95e0ddafb35` | `cmake/FindCHOLMOD.cmake` |
| `098ae524d0136556f7fd860df5cfc20e64d8e51feb4c7551ca345b5d69cb49c4` | `098ae524d0136556f7fd860df5cfc20e64d8e51feb4c7551ca345b5d69cb49c4` | `cmake/FindFFTW.cmake` |
| `2d1397cc46ff50700d06e62929af12a6c2454b10100d0463f75f6c063386efa0` | `2d1397cc46ff50700d06e62929af12a6c2454b10100d0463f75f6c063386efa0` | `cmake/FindLAPACKE.cmake` |
| `241b09e46cbda0bdc8c57c12d718c52ab8fd8148ac3d0d47cea7a9439273978b` | `241b09e46cbda0bdc8c57c12d718c52ab8fd8148ac3d0d47cea7a9439273978b` | `cmake/eigen-system/CMakeLists.txt` |
| `771ae4d540671ebe47c6bf97b78599bb8dd8a8eba8f0f08c581c66df1487f996` | `771ae4d540671ebe47c6bf97b78599bb8dd8a8eba8f0f08c581c66df1487f996` | `cmake/patches/stage02-fftwpp-odr.patch` |
| `86af407359b43072bf65e289ed92cdc20d677f4fd79c2c67767a03613757b33c` | `86af407359b43072bf65e289ed92cdc20d677f4fd79c2c67767a03613757b33c` | `cmake/patches/stage02-tomographymodels-odr.patch` |
| `3b9e1a584d1142f50e7fb0947ca82b0e6c55eba8ca7e30ecccd122c328d5fa2c` | `3b9e1a584d1142f50e7fb0947ca82b0e6c55eba8ca7e30ecccd122c328d5fa2c` | `examples/CMakeLists.txt` |
| `1adacfaca98b3ab5d62aceca2d51c1f2f28473fd9d1e5ceedf7c6b52a81b0dfd` | `1adacfaca98b3ab5d62aceca2d51c1f2f28473fd9d1e5ceedf7c6b52a81b0dfd` | `examples/OutputZiheng.cpp` |
| `ff69987bed879742c611ed6547340219efb9374c3fcb5ce01965f715d7e0daeb` | `ff69987bed879742c611ed6547340219efb9374c3fcb5ce01965f715d7e0daeb` | `examples/clean_bench_1.cpp` |
| `836cae39f74dfa1d6e64b66de4952f939fb18911d650568ea498635a94e5be58` | `836cae39f74dfa1d6e64b66de4952f939fb18911d650568ea498635a94e5be58` | `examples/clean_bench_10.cpp` |
| `5394cae84288372e5e693157d4fcb8236b8e19b85aaa1268426303aedf466504` | `5394cae84288372e5e693157d4fcb8236b8e19b85aaa1268426303aedf466504` | `examples/clean_bench_2.cpp` |
| `d78bd610c8515f5e0e51e8e095db90bc0c1d37e746520f4f595b12b981153cb4` | `d78bd610c8515f5e0e51e8e095db90bc0c1d37e746520f4f595b12b981153cb4` | `examples/clean_bench_3.cpp` |
| `827a76eeb48044ece5142a3ab77342e80cbdfbeec4589fdc85970c6a7aa4555a` | `827a76eeb48044ece5142a3ab77342e80cbdfbeec4589fdc85970c6a7aa4555a` | `examples/clean_bench_4.cpp` |
| `8df605a0df244ce52d1b2382a9dc7a3e57c978d5f52b9d0d8c6d8bded62c301d` | `8df605a0df244ce52d1b2382a9dc7a3e57c978d5f52b9d0d8c6d8bded62c301d` | `examples/clean_bench_5.cpp` |
| `c63b9ce360da4d29e42291213b9da04ec82b5efd5404d2f6261163413e1cf7a8` | `c63b9ce360da4d29e42291213b9da04ec82b5efd5404d2f6261163413e1cf7a8` | `examples/clean_bench_6.cpp` |
| `de052a8f4381f27037a0b2e99ad068593d8f6ac8f5d79ae57ad9f0334a852c89` | `de052a8f4381f27037a0b2e99ad068593d8f6ac8f5d79ae57ad9f0334a852c89` | `examples/clean_bench_7.cpp` |
| `5db1bfe83c509def987a79dc7b472d4be2a4a3af5acd94a85443b1e15dc85823` | `5db1bfe83c509def987a79dc7b472d4be2a4a3af5acd94a85443b1e15dc85823` | `examples/clean_bench_8.cpp` |
| `af66b7c6d0bbb1236376931298fd15e67b0558761950e6ecc84df9516ce88cfd` | `af66b7c6d0bbb1236376931298fd15e67b0558761950e6ecc84df9516ce88cfd` | `examples/clean_bench_9.cpp` |
| `f3a1deb714a2cc6f96d744ba8e245e60462e5dcaad37abe07db1c72e10af7eee` | `f3a1deb714a2cc6f96d744ba8e245e60462e5dcaad37abe07db1c72e10af7eee` | `examples/phobos_gravity.cpp` |
| `e81f985b2554107f432f07252d1b6796991a30789ecde505f33b6cdae7434694` | `e81f985b2554107f432f07252d1b6796991a30789ecde505f33b6cdae7434694` | `examples/phobos_gravity_perturbation.cpp` |
| `8f249fad456814f6dd416a124e5365dc35e1448cb6a727ee4128e04d3a51d45d` | `8f249fad456814f6dd416a124e5365dc35e1448cb6a727ee4128e04d3a51d45d` | `examples/phobos_heterogeneous.cpp` |
| `c699a78a3279812a2f0cb0c87a6274408914834bd3ba4caa763c5d6d69e056a0` | `c699a78a3279812a2f0cb0c87a6274408914834bd3ba4caa763c5d6d69e056a0` | `examples/wignertest.cpp` |
| `c8d6130345cf595d66ccb5ca2334503f1bab8d478d7d2469fc811ee93638712c` | `c8d6130345cf595d66ccb5ca2334503f1bab8d478d7d2469fc811ee93638712c` | `gplspec/All` |
| `60a010291acb8b5b6d2facd5595c1d92f7a35af4954b86aa8550587252ed6ddb` | `60a010291acb8b5b6d2facd5595c1d92f7a35af4954b86aa8550587252ed6ddb` | `gplspec/Core` |
| `59d5f04567e0f727cf26d2d532915f44e40f4b017085870b202cef22829a296c` | `59d5f04567e0f727cf26d2d532915f44e40f4b017085870b202cef22829a296c` | `gplspec/GeneralModels` |
| `4c108aba9c257dfddae459a8db6476edea54cbbecb5d6c11a05c2406ae8227af` | `4c108aba9c257dfddae459a8db6476edea54cbbecb5d6c11a05c2406ae8227af` | `gplspec/SimpleModels` |
| `99485412a4ac872226cd5e313102857068157fe0300dd20a0449a55d681cca07` | `99485412a4ac872226cd5e313102857068157fe0300dd20a0449a55d681cca07` | `gplspec/Test` |
| `5abf8304a5e8b2bc0c685d352c95e879560df592ef8fa79f4c428762f8081f28` | `5abf8304a5e8b2bc0c685d352c95e879560df592ef8fa79f4c428762f8081f28` | `gplspec/Timer` |
| `3a4532275d951fac0d83d7f8a5e7589013e4c7213c865de6f754e2d11c5e2975` | `3a4532275d951fac0d83d7f8a5e7589013e4c7213c865de6f754e2d11c5e2975` | `gplspec/src/GeneralModels/Density_Model_Constructors.h` |
| `838c211ee1539de51a3f56c4b5ff066a0df221afa5b091d152b9d42c2d719d68` | `838c211ee1539de51a3f56c4b5ff066a0df221afa5b091d152b9d42c2d719d68` | `gplspec/src/GeneralModels/Density_Model_Return.h` |
| `43ed635c487b79ca15ebacfa9bb57ca8b57a3e268d1b3b9460cf83aab088ffd3` | `43ed635c487b79ca15ebacfa9bb57ca8b57a3e268d1b3b9460cf83aab088ffd3` | `gplspec/src/GeneralModels/Earth_Density_Models_3D.h` |
| `b8cffa882f64af12308d01afafb8d2a805d0958a9188fba67cbb4b5fdfda6230` | `b8cffa882f64af12308d01afafb8d2a805d0958a9188fba67cbb4b5fdfda6230` | `gplspec/src/GeneralModels/Earth_General_Models_1D.h` |
| `71fd3884aa7e68f2de95578663cd43aed8492db7424e1e6ba1d2396c76001612` | `71fd3884aa7e68f2de95578663cd43aed8492db7424e1e6ba1d2396c76001612` | `gplspec/src/GeneralModels/Earth_General_Models_3D.h` |
| `c75f928843351c3199bbd549e43ffdef713e797891427819226e5e3e6592e935` | `c75f928843351c3199bbd549e43ffdef713e797891427819226e5e3e6592e935` | `gplspec/src/GeneralModels/Gravitational_Field.code-workspace` |
| `15b567a2d59727fec716e2890502c0b86d4835d3a8bd4b049dad56285458d5b6` | `15b567a2d59727fec716e2890502c0b86d4835d3a8bd4b049dad56285458d5b6` | `gplspec/src/GeneralModels/Gravity_Tools.h` |
| `7ec3167292f268f3b035cc9d22ecfe8819fbcca1a404f3e6ad1ab03bc57f7897` | `7ec3167292f268f3b035cc9d22ecfe8819fbcca1a404f3e6ad1ab03bc57f7897` | `gplspec/src/GeneralModels/Mapping_Perturbation_Utilities.h` |
| `13693634e7eb7d889bca05313158a11c99cfcd15ae3bd0332a45a493b863dc1c` | `13693634e7eb7d889bca05313158a11c99cfcd15ae3bd0332a45a493b863dc1c` | `gplspec/src/GeneralModels/Model_Construction_Utilities.h` |
| `d20f6604b96df1d7967acad2f502495ce3b32becbf0ea4c394240c1d2383d58c` | `d20f6604b96df1d7967acad2f502495ce3b32becbf0ea4c394240c1d2383d58c` | `gplspec/src/GeneralModels/Perturbation_Models.h` |
| `dd3ae893398648042c3cba03b6fe9a67d096cfdb283e474e8189c843d6073c06` | `dd3ae893398648042c3cba03b6fe9a67d096cfdb283e474e8189c843d6073c06` | `gplspec/src/Gravitational_Field.code-workspace` |
| `e7dbcc0735894e4231f343c18c30aeeb0ca1078d7fa7ca469730bd28cbee0547` | `e7dbcc0735894e4231f343c18c30aeeb0ca1078d7fa7ca469730bd28cbee0547` | `gplspec/src/Mapping_Tools.h` |
| `7ca5f3acf5a5b56f1952f750fecc253a731afa99231e3cdc801de609841f0f8b` | `7ca5f3acf5a5b56f1952f750fecc253a731afa99231e3cdc801de609841f0f8b` | `gplspec/src/Pseudospectral_Matrix_Wrapper.h` |
| `21b949b669e118c33a1e8c80d3e45ef99f2e7418839d9e2a723809d832ab1123` | `21b949b669e118c33a1e8c80d3e45ef99f2e7418839d9e2a723809d832ab1123` | `gplspec/src/Pseudospectral_Matrix_Wrapper3D.h` |
| `8ecb567e9bc0a8d295c912c7b1db78e84573fccb921c6b181f3b1105d2a8a8f5` | `8ecb567e9bc0a8d295c912c7b1db78e84573fccb921c6b181f3b1105d2a8a8f5` | `gplspec/src/Radial_Tools.h` |
| `ab34188526a0d495fd8bbe4315417f17d42b1fc619ed19b169f39e61f6761e35` | `ab34188526a0d495fd8bbe4315417f17d42b1fc619ed19b169f39e61f6761e35` | `gplspec/src/Shared_Utilities.h` |
| `16f44a6ebf29ff3010700fbb798fb00fd33b3d8874c27b1825c158c0fe489e34` | `16f44a6ebf29ff3010700fbb798fb00fd33b3d8874c27b1825c158c0fe489e34` | `gplspec/src/SimpleModels/Simple_Model_Class_Definition.h` |
| `f107f13163eb4b6511562dd0ba7c50f1f6c84b4252aa89672eef3a542ed76354` | `f107f13163eb4b6511562dd0ba7c50f1f6c84b4252aa89672eef3a542ed76354` | `gplspec/src/SimpleModels/Simple_Model_Constructor_Definitions.h` |
| `8db29030e8b862f507142959b5310006ce884528da240beb851786364746bfba` | `8db29030e8b862f507142959b5310006ce884528da240beb851786364746bfba` | `gplspec/src/SimpleModels/Simple_Model_Return_Definitions.h` |
| `901fe7bc431ecab8c8bc1a8935459daccfcc6301f65065fbca5ec62a1598d75f` | `901fe7bc431ecab8c8bc1a8935459daccfcc6301f65065fbca5ec62a1598d75f` | `gplspec/src/Spectral_Element_Tools.h` |
| `b51e29135a2b2a33dfbec6aa11ba96f48163b978a79fe3f753724da745ed31f2` | `b51e29135a2b2a33dfbec6aa11ba96f48163b978a79fe3f753724da745ed31f2` | `gplspec/src/SphericalGeometryPreconditioner.h` |
| `a3c380a963301c7eba27edc9d65ca065e3123d305a25ec0433463a959c225450` | `a3c380a963301c7eba27edc9d65ca065e3123d305a25ec0433463a959c225450` | `gplspec/src/Spherical_Integrator.h` |
| `5433fcbc5ee0124797266db07649afc45039b75e9cf8d03e1821fb7be50e4c93` | `5433fcbc5ee0124797266db07649afc45039b75e9cf8d03e1821fb7be50e4c93` | `gplspec/src/Timer_Class.h` |
| `c9f8532616e25c19524e87691d3a71dcffcc9955fa9984dadfd5ddaaf631f995` | `c9f8532616e25c19524e87691d3a71dcffcc9955fa9984dadfd5ddaaf631f995` | `gplspec/src/testremove/Aspherical_Geometry_Potential.cpp` |
| `8d52e89601786b35555bdfe01f80086e98baa6d81404d62c56fea043a3246b81` | `8d52e89601786b35555bdfe01f80086e98baa6d81404d62c56fea043a3246b81` | `gplspec/src/testremove/Spectral_Element.cpp` |
| `6743e5bf3b44465f7a407d8c2ea271324ec5333b45b5e82dc30a708eedc2f73a` | `6743e5bf3b44465f7a407d8c2ea271324ec5333b45b5e82dc30a708eedc2f73a` | `gplspec/src/testremove/Spectral_Element_Sparse.cpp` |
| `fc109362afb3a78e8dca034977c4abb98814508ce0504b5781fc852fa6b18a4d` | `fc109362afb3a78e8dca034977c4abb98814508ce0504b5781fc852fa6b18a4d` | `gplspec/src/testremove/Spectral_Element_Sparse_Test.cpp` |
| `a2d12cb9053212381e8d0644305a00cc411f2ff48fdf62d937cf01638b3b9926` | `a2d12cb9053212381e8d0644305a00cc411f2ff48fdf62d937cf01638b3b9926` | `gplspec/src/testremove/Spherical_Geometry_Potential.cpp` |
| `4cc07945083ce8aeb7455375c7e5bfbf9d78bb45d9370e91ab27f234112a390f` | `4cc07945083ce8aeb7455375c7e5bfbf9d78bb45d9370e91ab27f234112a390f` | `gplspec/src/testremove/Test.cpp` |
| `a79a61aa3f30d89e89ed49ab8e8b863a471fe9f4eec4bcc183c66b99782a3763` | `a79a61aa3f30d89e89ed49ab8e8b863a471fe9f4eec4bcc183c66b99782a3763` | `gplspec/src/testtools.h` |
| `b7011661ac90d71c7f655b4ffa8d4e1dfa83f458493759c2f471f63f0ac52b34` | `350c39aab2e4407bc886405f2e4ee495d83947e68029851c27a029091a701c1e` | `tests/CMakeLists.txt` |
| `036c64845937c5ee87f6aa91af3809d8df10c92c9797ae8a16132abcb7202268` | `036c64845937c5ee87f6aa91af3809d8df10c92c9797ae8a16132abcb7202268` | `tests/compare_reference.py` |
| `UNTRACKED` | `1152e096d24060e9918f4dc4ff45bfa64e8fe1e5e35b60556ec36de776e81c18` | `tests/external_consumer/CMakeLists.txt` |
| `UNTRACKED` | `93ff5d621390d8f78149ccb57001e7d6ee07967fc84589b479ea0cf336dee852` | `tests/external_consumer/main.cpp` |
| `UNTRACKED` | `7c227c1cb40a4dc28a9b3b4cb83477e9d72ca11464fcf66cf07e0c47fe8d0aba` | `tests/external_consumer/peer.cpp` |
| `8a0949744ec1a89918e653be38f82dae07b5aa93d9e2d837474c509b1845b809` | `8a0949744ec1a89918e653be38f82dae07b5aa93d9e2d837474c509b1845b809` | `tests/extract_tutorial_cpp.py` |
| `42028cbc11754dde6593e4e3e25346c659414b209d5849179553ef8d83f785e7` | `42028cbc11754dde6593e4e3e25346c659414b209d5849179553ef8d83f785e7` | `tests/pr07_rotated_output.cpp` |
| `UNTRACKED` | `103d812a30d720ad15e7af5175203fc89faf54053c355baf450839cd3b9681c0` | `tests/run_external_consumer.py` |
| `9afaec3ef0188c3d233cb249b9ab9214218a75ee111d0d50019c0be6ea64f5de` | `9afaec3ef0188c3d233cb249b9ab9214218a75ee111d0d50019c0be6ea64f5de` | `tests/run_stage00.sh` |
| `2c133aa19b3e67a1c14d24e1072555cde56de5769b4162c37e9d121354c48c54` | `2c133aa19b3e67a1c14d24e1072555cde56de5769b4162c37e9d121354c48c54` | `tests/run_stage05_1d_reference.sh` |
| `e9adeeed9b9c5dbec5ddccbf37a1014a3ee1414b271734295751fd987b403270` | `e9adeeed9b9c5dbec5ddccbf37a1014a3ee1414b271734295751fd987b403270` | `tests/run_stage05_file_perturbation_reference.sh` |
| `babf13a6b7c7ba7eed2a9052626d42a430df47fbb14c6a1b3565ef8a53b2c775` | `babf13a6b7c7ba7eed2a9052626d42a430df47fbb14c6a1b3565ef8a53b2c775` | `tests/run_stage06_output_reference.py` |
| `4f3d141a837da4efc6605fd99c9c99d7c4f956829fc7108dcaf5f32e23c59e78` | `4f3d141a837da4efc6605fd99c9c99d7c4f956829fc7108dcaf5f32e23c59e78` | `tests/stage00_reference.cpp` |
| `27aff5ddd8eb1cbae16d5ac345a3f45b364965501b6661a30c71393e2e6dee79` | `27aff5ddd8eb1cbae16d5ac345a3f45b364965501b6661a30c71393e2e6dee79` | `tests/stage02_all_smoke.cpp` |
| `737c7c4a706553d470da33fcd6415a57723e6b8ceb4899680d8919e601bc3a03` | `737c7c4a706553d470da33fcd6415a57723e6b8ceb4899680d8919e601bc3a03` | `tests/stage02_core_smoke.cpp` |
| `f79413eaee9070d41ffccca1e8ca5d459e2818e7a1624bf865edcfea7aa8d866` | `f79413eaee9070d41ffccca1e8ca5d459e2818e7a1624bf865edcfea7aa8d866` | `tests/stage02_generalmodels_smoke.cpp` |
| `3233eae59aecbb6bf506fb32b37fc47f3ab6900330a0cef61783dbe42c1c1c29` | `3233eae59aecbb6bf506fb32b37fc47f3ab6900330a0cef61783dbe42c1c1c29` | `tests/stage02_header_main.cpp` |
| `605f9386ea9b79e84e27285eed032559736b57171d08a31d25d49d778b8c1f4c` | `605f9386ea9b79e84e27285eed032559736b57171d08a31d25d49d778b8c1f4c` | `tests/stage02_header_peer.cpp` |
| `5e4b62699836c129abceea69ffa5d97b5ada689f6f6168889beb77d3a5694fbb` | `5e4b62699836c129abceea69ffa5d97b5ada689f6f6168889beb77d3a5694fbb` | `tests/stage02_simplemodels_smoke.cpp` |
| `UNTRACKED` | `25e972cd622e0a3bef37f83e341966b97e17e0cb8eaeede6962d70af48cd4439` | `tests/stage02_test_smoke.cpp` |
| `UNTRACKED` | `aa04c4cea1bc41506aac9e2fb92c9beb2121414411cd4b8e29b9f161a51993b7` | `tests/stage02_timer_smoke.cpp` |
| `7651fbfff9ac14599086faaad6b0bfa1401449f466cfdeb975a27b946e384e14` | `7651fbfff9ac14599086faaad6b0bfa1401449f466cfdeb975a27b946e384e14` | `tests/stage03_utilities.cpp` |
| `1afea93feff2bafd6392719b109a43a886468c8a3927af87c63782948a4188eb` | `1afea93feff2bafd6392719b109a43a886468c8a3927af87c63782948a4188eb` | `tests/stage04_construction_storage.cpp` |
| `41013766a4ff87b9f0d810bd2c696caa64171dcf04abd2dcc213d7bdd1432702` | `41013766a4ff87b9f0d810bd2c696caa64171dcf04abd2dcc213d7bdd1432702` | `tests/stage04_density_sampling.cpp` |
| `3052df1403c21db60254b0cbcdcd47c5c2d3fa1714352e641346958c6e743126` | `3052df1403c21db60254b0cbcdcd47c5c2d3fa1714352e641346958c6e743126` | `tests/stage04_discretization_setup.cpp` |
| `d3b7473fcb96068ca3e12a4e6aa797d203522d3b1a4173d2ae2e7e8a27436411` | `d3b7473fcb96068ca3e12a4e6aa797d203522d3b1a4173d2ae2e7e8a27436411` | `tests/stage04_mapping_geometry.cpp` |
| `5103f624ad62f8ceb3317668c1d1caf21a248cacb7160a104c209cb9c6f5740c` | `5103f624ad62f8ceb3317668c1d1caf21a248cacb7160a104c209cb9c6f5740c` | `tests/stage05_file_perturbation_reference.cpp` |
| `15792d9630fbc8a328d3bc89f6c706ed0f566be695f2f12f067d332ee94c48ad` | `15792d9630fbc8a328d3bc89f6c706ed0f566be695f2f12f067d332ee94c48ad` | `tests/stage05_operator_reference.cpp` |
| `b2a4bcc70a8c724e49c0fed3612280b0f703d1dc067096ea9440623e45238f37` | `b2a4bcc70a8c724e49c0fed3612280b0f703d1dc067096ea9440623e45238f37` | `tests/stage05_tensor_contraction.cpp` |
| `9b44e5d2a7e5cd2be7ede39676a58a3d491b52981bf7ffde700fc20fd1eac5ba` | `9b44e5d2a7e5cd2be7ede39676a58a3d491b52981bf7ffde700fc20fd1eac5ba` | `tests/stage06_output_helpers.cpp` |
| `843e8f70eec0f34c3d86da7d5064ca3ff4375979d6faa39b4960302efafa0425` | `843e8f70eec0f34c3d86da7d5064ca3ff4375979d6faa39b4960302efafa0425` | `tests/stage06_output_reference.cpp` |
| `UNTRACKED` | `87d958f306c9ebae4540eb8ccea337713f7164a58f9affab7da3a3192236bc2a` | `tests/stage1c_file_tomography_solver.cpp` |
| `33d2787e6104c05492a5c6adaec1a07018dc3167219b3ac9f2910f7c4e99b42b` | `33d2787e6104c05492a5c6adaec1a07018dc3167219b3ac9f2910f7c4e99b42b` | `tests/test_stage02_patch_application.py` |

## Inputs and frozen references

| SHA-256 | Path |
| --- | --- |
| `ec9869508847eb6753cf0165f6176dd3eea1f03ac5717533ffb6dcfe09489685` | `modeldata/PREM500` |
| `4ebd0ca09731a593c7960ced6a8ddde0abbe03fa91e49e70cf897dd41019c73e` | `modeldata/S40RTS_dvs.nc` |
| `a902254d3d416c55f103058c045a30036e7b5edd8b56166e086d43de4de13c3f` | `tests/reference/homogeneous-MatrixSolution.out` |
| `ef20bd7e3016661c60903c290c74c599151fdfc0417443bcf81c4215401f105e` | `tests/reference/stage00-baseline.csv` |
| `a9e89b9ac1250f3f01b695ec5be0fcb79990a85da42e6eb80ae001f5f88a43d6` | `tests/reference/stage05-1d-original.csv` |
| `595527f990abdbc2e32e31cf289915ec3c2cc647733609f86462282d68dca3e9` | `tests/reference/stage05-file-perturbation-coefficients.dat` |
| `08cf04f1243fa128645490c4eba13eed5f7b5fceff76b032f2f5d8141efe3263` | `tests/reference/stage05-file-perturbation-original.csv` |
| `f58a90f7f05c1f6e88712a8d6b20f4202284cb2bb266782961e1a018c8b22236` | `tests/reference/stage06-output-baseline.sha256` |
