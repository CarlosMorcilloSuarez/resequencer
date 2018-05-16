#!/bin/bash

# FUNCTIONS ===============================================================

# Executes the given command and checks if the output file is
# created with the right name and if its contents are identical
# to the model file
test_created_files(){
  TEST_NAME=$1
  COMMAND=$2
  TARGET_DIRECTORY_NAME=$3
  REFERENCE_DIRECTORY_NAME=$4

  echo
  echo
  echo ${TEST_NAME}
  echo

  # Cleans target directory before command execution
  rm -r ./testingFiles/${TARGET_DIRECTORY_NAME}/*

  # Executes command
  eval ${COMMAND}

  # Checks whether both directories have identical content
  targetContent=$( \
      find ./testingFiles/${TARGET_DIRECTORY_NAME} -type f \
      | sed -r s/^\\.\\/testingFiles\\/${TARGET_DIRECTORY_NAME}//g \
      | sed s/^\\///g )
  referenceContent=$( \
      find ./testingFiles/${REFERENCE_DIRECTORY_NAME} -type f \
      | sed -r s/^\\.\\/testingFiles\\/${REFERENCE_DIRECTORY_NAME}//g \
      | sed s/^\\///g )

  STATUS='OK'
  for file in $(echo "${targetContent} ${referenceContent}" \
                      | tr ' ' '\n' \
                      | sort -u)
  do
  	ZDIFF_OUTPUT=$(zdiff ./testingFiles/${TARGET_DIRECTORY_NAME}/${file} \
                         ./testingFiles/${REFERENCE_DIRECTORY_NAME}/${file})
    if test $? -ne 0
    then
      STATUS='NOK'
      echo ./testingFiles/${TARGET_DIRECTORY_NAME}/${file}
      echo "${ZDIFF_OUTPUT}"
    else
      echo ./testingFiles/${TARGET_DIRECTORY_NAME}/${file}
      echo "OK"
    fi
  done

  if test $STATUS = 'OK'
  then
    echo "------------------------------------------------- OK"
  else
    echo "====================================================== FAIL"
  fi
}




# TESTS ====================================================================

# Test -----------------------------------------------------------------
# Config Area
TEST_NAME="test_fastResequencer_default_values"
COMMAND='
      python fastResequencer.py
      --reference ./testingFiles/Rat-Monkey_Reference.fa
      --output-file ./testingFiles/fastTest1/fastTest1.fq
      --seed 1969
      '
# The previous command should create output files in the directory:
TARGET_DIRECTORY_NAME='fastTest1'
# That should be identical to those in the directory:
REFERENCE_DIRECTORY_NAME='fastTest1_Reference'


# Execution Area - don't touch it
test_created_files "${TEST_NAME}" \
                  "${COMMAND}" \
                  "${TARGET_DIRECTORY_NAME}" \
                  "${REFERENCE_DIRECTORY_NAME}"



# Test -----------------------------------------------------------------
# Config Area
TEST_NAME="test_fastResequencer_default_values_zip"
COMMAND='
      python fastResequencer.py
      --reference ./testingFiles/Rat-Monkey_Reference.fa
      --output-file ./testingFiles/fastTest2/sample.fq.gz
      --seed 1969
      '
# The previous command should create output files in the directory:
TARGET_DIRECTORY_NAME='fastTest2'
# That should be identical to those in the directory:
REFERENCE_DIRECTORY_NAME='fastTest2_Reference'


# Execution Area - don't touch it
test_created_files "${TEST_NAME}" \
                  "${COMMAND}" \
                  "${TARGET_DIRECTORY_NAME}" \
                  "${REFERENCE_DIRECTORY_NAME}"
