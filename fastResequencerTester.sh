#!/bin/bash


run_tests(){

# Tests are written here ###################################

# Test ----------------------------------------
initialize_test "Generates fastq file with default_values"

python fastResequencer.py \
  --reference ./test/inputs/Rat-Monkey_Reference.fa \
  --output-file ./test/outputs/fastTest1.fq \
  --seed 1969

compare_files \
  ./test/outputs/fastTest1.fq \
  ./test/inputs/references/Model_fastTest1.fq


# Test ----------------------------------------
initialize_test "Generates gzip fastq file"

python fastResequencer.py \
  --reference ./test/inputs/Rat-Monkey_Reference.fa \
  --output-file ./test/outputs/fastTest2.fq.gz \
  --seed 1969

compare_files \
  ./test/outputs/fastTest2.fq.gz \
  ./test/inputs/references/Model_fastTest2.fq.gz


  # Test ----------------------------------------
  initialize_test "Fragments Genome into 100 length, 25 overlap reads"

  python fastResequencer.py \
    --reference ./test/inputs/Rat-Monkey_Reference.fa \
    --output-file ./test/outputs/fastTest3.fq \
    --full 25 \
    --length 100 \
    --seed 1969

  compare_files \
    ./test/outputs/fastTest3.fq \
    ./test/inputs/references/Model_fastTest3.fq

############################################################
}


# Testing Functions

initialize_test(){
  echo
  echo $1
}


compare_strings(){
  STRING=$1
  REFERENCE_STRING=$2

  if [ "${STRING}" = "${REFERENCE_STRING}" ]
  then
    echo "------------------------------------------------- OK"
  else
    echo "ERROR ------------------------------ ERROR"
    echo ${STRING}---${REFERENCE_STRING}--- Seem to be different
  fi
}


compare_files(){
  FILE_NAME=$1
  MODEL_FILE=$2

  # Checks expected outputfile
  if [ -f "${FILE_NAME}" ]
  then
      zdiff ${FILE_NAME} ${MODEL_FILE} > /dev/null
      if test $? -eq 0
      then
        echo "------------------------------------------------- OK"
      else
        echo "ERROR ------------------------------ ERROR"
        echo ${FILE_NAME} --- ${MODEL_FILE}     --- Seem to be different
      fi
  else
      echo 'ERROR ------------------------------ ERROR'
      echo "Expected ${FILE_NAME} not found."
  fi
}



compare_directories(){
  TARGET_DIRECTORY_NAME=$1
  REFERENCE_DIRECTORY_NAME=$2

  # Checks whether both directories have identical content
  TARGET_DIRECTORY_NAME_LENGTH=${#TARGET_DIRECTORY_NAME}
  targetContent=$(
                find ${TARGET_DIRECTORY_NAME} -type f |
                cut -c$((${TARGET_DIRECTORY_NAME_LENGTH}+2))-
                )

  REFERENCE_DIRECTORY_NAME_LENGTH=${#REFERENCE_DIRECTORY_NAME}
  referenceContent=$(
                find ${REFERENCE_DIRECTORY_NAME} -type f |
                cut -c$((${REFERENCE_DIRECTORY_NAME_LENGTH}+2))-
                )

  STATUS='OK'
  for file in $(echo "${targetContent} ${referenceContent}" \
                      | tr ' ' '\n' \
                      | sort -u)
  do
    ZDIFF_OUTPUT=$(zdiff ${TARGET_DIRECTORY_NAME}/${file} \
                        ${REFERENCE_DIRECTORY_NAME}/${file})
    if test $? -ne 0
    then
      STATUS='NOK'
      echo -e '\t'${TARGET_DIRECTORY_NAME}/${file}
      echo -e '\t'"${ZDIFF_OUTPUT}"
    #else
      #echo -e '\t'${TARGET_DIRECTORY_NAME}/${file}
      #echo -e '\t'"OK"
    fi
  done

  if test $STATUS = 'OK'
  then
    echo "------------------------------------------------- OK"
  else
    echo "====================================================== FAIL"
  fi
}



# Creates and cleans test directory structure
if true; then
  if [ ! -d "test" ]; then

    mkdir test
    mkdir test/scripts
    mkdir test/inputs
    mkdir test/tmp
    mkdir test/outputs

  else

    if [ ! -d "test/scripts" ]; then
    mkdir test/scripts
    fi

    if [ ! -d "test/inputs" ]; then
    mkdir test/inputs
    fi

    if [ ! -d "test/tmp" ]; then
    mkdir test/tmp
    else
    rm -fR test/tmp/*
    fi

    if [ ! -d "test/outputs" ]; then
    mkdir test/outputs
    else
    rm -fR test/outputs/*
    fi

  fi
else
  echo "WARNING -- TESTING MODE"
fi

run_tests
