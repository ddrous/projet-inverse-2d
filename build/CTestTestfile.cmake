# CMake generated Testfile for 
# Source directory: /mnt/c/Users/Roussel/Dropbox/Unistra/SEMESTRE 2/Projet & Stage/Inverse/REPO_2D
# Build directory: /mnt/c/Users/Roussel/Dropbox/Unistra/SEMESTRE 2/Projet & Stage/Inverse/REPO_2D/build
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(simple_run "/mnt/c/Users/Roussel/Dropbox/Unistra/SEMESTRE 2/Projet & Stage/Inverse/REPO_2D/build/transfer" "../src/config/test_1.cfg")
add_test(no_config "/mnt/c/Users/Roussel/Dropbox/Unistra/SEMESTRE 2/Projet & Stage/Inverse/REPO_2D/build/transfer")
set_tests_properties(no_config PROPERTIES  WILL_FAIL "TRUE")
add_test(bad_param "/mnt/c/Users/Roussel/Dropbox/Unistra/SEMESTRE 2/Projet & Stage/Inverse/REPO_2D/build/transfer" "../src/config/test_2.cfg")
set_tests_properties(bad_param PROPERTIES  WILL_FAIL "TRUE")
add_test(unparsable "/mnt/c/Users/Roussel/Dropbox/Unistra/SEMESTRE 2/Projet & Stage/Inverse/REPO_2D/build/transfer" "../src/config/test_3.cfg")
set_tests_properties(unparsable PROPERTIES  WILL_FAIL "TRUE")
