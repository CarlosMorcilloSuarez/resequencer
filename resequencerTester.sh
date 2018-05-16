#!/bin/bash

# FUNCTIONS ===============================================================

# Executes the given command and checks if the output file is
# created with the right name and if its contents are identical
# to the model file
test_file_created(){
  TEST_NAME=$1
  COMMAND=$2
  FILE_NAME=$3
  MODEL_FILE=$4

  echo
  echo ${TEST_NAME}

  # Executes command
  ${COMMAND}

  # Checks expected outputfile
  if [ -f "${FILE_NAME}" ]
  then
      diff ${FILE_NAME} ${MODEL_FILE} > /dev/null
      if test $? -eq 0
      then
        echo "OK"
      else
        echo "ERROR ------------------------------ ERROR"
        echo ${FILE_NAME} --- ${MODEL_FILE}     --- Seem to be different
      fi
  else
      echo 'ERROR ------------------------------ ERROR'
	    echo "Expected ${FILE_NAME} not found."
  fi

  # Removes output file
  if [ -f "${FILE_NAME}" ]
  then
      rm ${FILE_NAME}
  fi

}


# TESTS ====================================================================

# Test -----------------------------------------------------------------------
# Config Area
TEST_NAME="test_resequencer_default_values"
COMMAND='
      python resequencer.py
      --reference ./testingFiles/Rat-Monkey_Reference.fa
      --name ./testingFiles/Test1
      --seed 1969
      '
# The previous command should generate an output file with name:
FILE_NAME='./testingFiles/Test1.fastq'
# Whose contents should be identical to the file:
MODEL_FILE='./testingFiles/Model_Test1.fastq'


# Execution Area - don't touch it
test_file_created "${TEST_NAME}" "${COMMAND}" "${FILE_NAME}" "${MODEL_FILE}"




# Test -----------------------------------------------------------------------
# Config Area
TEST_NAME="test_resequencer_coverage_error-rate"
COMMAND='
      python resequencer.py
            --reference ./testingFiles/Rat-Monkey_Reference.fa
            --name ./testingFiles/Test2
            --coverage 5
            --error-rate 0.1
            --seed 1969
      '
# The previous command should generate an output file with name:
FILE_NAME='./testingFiles/Test2.fastq'
# Whose contents should be identical to the file:
MODEL_FILE='./testingFiles/Model_Test2.fastq'


# Execution Area - don't touch it
test_file_created "${TEST_NAME}" "${COMMAND}" "${FILE_NAME}" "${MODEL_FILE}"


# Test -----------------------------------------------------------------------
# Config Area
TEST_NAME="test_resequencer_cnv"
COMMAND='
      python resequencer.py
            --reference ./testingFiles/Rat-Monkey_Reference.fa
            --name ./testingFiles/Test3
            --seed 1969
            --duplications ./testingFiles/CNV3.conf
      '
# The previous command should generate an output file with name:
FILE_NAME='./testingFiles/Test3.fastq'
# Whose contents should be identical to the file:
MODEL_FILE='./testingFiles/Model_Test3.fastq'


# Execution Area - don't touch it
test_file_created "${TEST_NAME}" "${COMMAND}" "${FILE_NAME}" "${MODEL_FILE}"


# Test -----------------------------------------------------------------------
# Config Area
TEST_NAME="test_resequencer_read-length"
COMMAND='
      python resequencer.py
            --reference ./testingFiles/Rat-Monkey_Reference.fa
            --name ./testingFiles/Test4
            --seed 1969
            --length 25
      '
# The previous command should generate an output file with name:
FILE_NAME='./testingFiles/Test4.fastq'
# Whose contents should be identical to the file:
MODEL_FILE='./testingFiles/Model_Test4.fastq'


# Execution Area - don't touch it
test_file_created "${TEST_NAME}" "${COMMAND}" "${FILE_NAME}" "${MODEL_FILE}"
