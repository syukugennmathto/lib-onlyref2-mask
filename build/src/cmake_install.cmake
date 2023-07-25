# Install script for directory: /Users/s-yokoyama_lab/Downloads/liboqs-main-2/src

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/Library/Developer/CommandLineTools/usr/bin/objdump")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/liboqs" TYPE FILE FILES
    "/Users/s-yokoyama_lab/Downloads/liboqs-main-2/build/src/liboqsConfig.cmake"
    "/Users/s-yokoyama_lab/Downloads/liboqs-main-2/build/src/liboqsConfigVersion.cmake"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/pkgconfig" TYPE FILE FILES "/Users/s-yokoyama_lab/Downloads/liboqs-main-2/build/src/liboqs.pc")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/Users/s-yokoyama_lab/Downloads/liboqs-main-2/build/lib/liboqs.a")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/liboqs.a" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/liboqs.a")
    execute_process(COMMAND "/Library/Developer/CommandLineTools/usr/bin/ranlib" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/liboqs.a")
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/liboqs/liboqsTargets.cmake")
    file(DIFFERENT _cmake_export_file_changed FILES
         "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/liboqs/liboqsTargets.cmake"
         "/Users/s-yokoyama_lab/Downloads/liboqs-main-2/build/src/CMakeFiles/Export/c7e97583fbc7c9ca02085e7795e05761/liboqsTargets.cmake")
    if(_cmake_export_file_changed)
      file(GLOB _cmake_old_config_files "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/liboqs/liboqsTargets-*.cmake")
      if(_cmake_old_config_files)
        string(REPLACE ";" ", " _cmake_old_config_files_text "${_cmake_old_config_files}")
        message(STATUS "Old export file \"$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/liboqs/liboqsTargets.cmake\" will be replaced.  Removing files [${_cmake_old_config_files_text}].")
        unset(_cmake_old_config_files_text)
        file(REMOVE ${_cmake_old_config_files})
      endif()
      unset(_cmake_old_config_files)
    endif()
    unset(_cmake_export_file_changed)
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/liboqs" TYPE FILE FILES "/Users/s-yokoyama_lab/Downloads/liboqs-main-2/build/src/CMakeFiles/Export/c7e97583fbc7c9ca02085e7795e05761/liboqsTargets.cmake")
  if(CMAKE_INSTALL_CONFIG_NAME MATCHES "^()$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/liboqs" TYPE FILE FILES "/Users/s-yokoyama_lab/Downloads/liboqs-main-2/build/src/CMakeFiles/Export/c7e97583fbc7c9ca02085e7795e05761/liboqsTargets-noconfig.cmake")
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/oqs" TYPE FILE FILES
    "/Users/s-yokoyama_lab/Downloads/liboqs-main-2/src/oqs.h"
    "/Users/s-yokoyama_lab/Downloads/liboqs-main-2/src/common/common.h"
    "/Users/s-yokoyama_lab/Downloads/liboqs-main-2/src/common/rand/rand.h"
    "/Users/s-yokoyama_lab/Downloads/liboqs-main-2/src/common/aes/aes.h"
    "/Users/s-yokoyama_lab/Downloads/liboqs-main-2/src/common/sha2/sha2.h"
    "/Users/s-yokoyama_lab/Downloads/liboqs-main-2/src/common/sha3/sha3.h"
    "/Users/s-yokoyama_lab/Downloads/liboqs-main-2/src/common/sha3/sha3x4.h"
    "/Users/s-yokoyama_lab/Downloads/liboqs-main-2/src/kem/kem.h"
    "/Users/s-yokoyama_lab/Downloads/liboqs-main-2/src/sig/sig.h"
    "/Users/s-yokoyama_lab/Downloads/liboqs-main-2/src/kem/bike/kem_bike.h"
    "/Users/s-yokoyama_lab/Downloads/liboqs-main-2/src/kem/frodokem/kem_frodokem.h"
    "/Users/s-yokoyama_lab/Downloads/liboqs-main-2/src/kem/ntruprime/kem_ntruprime.h"
    "/Users/s-yokoyama_lab/Downloads/liboqs-main-2/src/kem/classic_mceliece/kem_classic_mceliece.h"
    "/Users/s-yokoyama_lab/Downloads/liboqs-main-2/src/kem/hqc/kem_hqc.h"
    "/Users/s-yokoyama_lab/Downloads/liboqs-main-2/src/kem/kyber/kem_kyber.h"
    "/Users/s-yokoyama_lab/Downloads/liboqs-main-2/src/sig/dilithium/sig_dilithium.h"
    "/Users/s-yokoyama_lab/Downloads/liboqs-main-2/src/sig/falcon/sig_falcon.h"
    "/Users/s-yokoyama_lab/Downloads/liboqs-main-2/src/sig/sphincs/sig_sphincs.h"
    "/Users/s-yokoyama_lab/Downloads/liboqs-main-2/build/include/oqs/oqsconfig.h"
    )
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/Users/s-yokoyama_lab/Downloads/liboqs-main-2/build/src/common/cmake_install.cmake")
  include("/Users/s-yokoyama_lab/Downloads/liboqs-main-2/build/src/kem/bike/cmake_install.cmake")
  include("/Users/s-yokoyama_lab/Downloads/liboqs-main-2/build/src/kem/frodokem/cmake_install.cmake")
  include("/Users/s-yokoyama_lab/Downloads/liboqs-main-2/build/src/kem/ntruprime/cmake_install.cmake")
  include("/Users/s-yokoyama_lab/Downloads/liboqs-main-2/build/src/kem/classic_mceliece/cmake_install.cmake")
  include("/Users/s-yokoyama_lab/Downloads/liboqs-main-2/build/src/kem/hqc/cmake_install.cmake")
  include("/Users/s-yokoyama_lab/Downloads/liboqs-main-2/build/src/kem/kyber/cmake_install.cmake")
  include("/Users/s-yokoyama_lab/Downloads/liboqs-main-2/build/src/sig/dilithium/cmake_install.cmake")
  include("/Users/s-yokoyama_lab/Downloads/liboqs-main-2/build/src/sig/falcon/cmake_install.cmake")
  include("/Users/s-yokoyama_lab/Downloads/liboqs-main-2/build/src/sig/sphincs/cmake_install.cmake")

endif()

