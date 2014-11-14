# Copyright (C) 2008-2014 LAAS-CNRS, JRL AIST-CNRS.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# DISTCHECK_SETUP
# ---------------
#
# Add a distcheck target to check the generated tarball.
#
# This step calls `make distdir' to generate a copy of the project without
# the git history and with the `.version' file (as it will be when an user
# will retrieve a stable version).
# Then:
# - create _build and _inst to respectively create a build and an installation
#   directory.
# - copy the CMakeCache.txt file.
# - run cmake with _inst as the installation prefix
# - run make, make check, make install and make uninstall
# - remove _build and _inst.
#
# During the compilation phase, all files in the source tree are modified
# to *not* be writeable to detect bad compilation steps which tries to modify
# the source tree. Permissions are reverted at the end of the check.
#
MACRO(DISTCHECK_SETUP)
  FIND_PROGRAM(SED sed)
  FIND_PROGRAM(TAR tar)
  FIND_PROGRAM(GZIP gzip)
  SET(INSTDIR ${CMAKE_BINARY_DIR}/${PACKAGE_NAME}-${PACKAGE_VERSION}/_inst)

  ADD_CUSTOM_TARGET(dist 
    COMMAND
    cd ${CMAKE_SOURCE_DIR}
    && git archive --format=tar --prefix=${PACKAGE_NAME}-${PACKAGE_VERSION}/ HEAD | 
    ${GZIP} > ${CMAKE_BINARY_DIR}/${PACKAGE_NAME}-${PACKAGE_VERSION}.tar.gz
  )

#add_custom_target(distcheck cd ${CMAKE_BINARY_DIR} &&
#                        rm -rf ${PACKAGE_NAME}-${PACKAGE_VERSION} &&
#                        gzip -d ${PACKAGE_NAME}-${PACKAGE_VERSION}.tar.gz &&
#                        tar -xf ${PACKAGE_NAME}-${PACKAGE_VERSION}.tar &&
#                        cd ${PACKAGE_NAME}-${PACKAGE_VERSION}/ &&
#                        mkdir _build && cd _build &&
#                        cmake .. && make && make test &&
#                        cd ${CMAKE_BINARY_DIR} &&
#                        tar -rf ${PACKAGE_NAME}-${PACKAGE_VERSION}.tar ${PACKAGE_NAME}-${PACKAGE_VERSION}/doc/ &&
#                        gzip ${CMAKE_BINARY_DIR}/${PACKAGE_NAME}-${PACKAGE_VERSION}.tar)
#add_dependency(distcheck dist)

#    cd ${CMAKE_BINARY_DIR}
  ADD_CUSTOM_TARGET(distcheck
    COMMAND
    rm -rf ${PACKAGE_NAME}-${PACKAGE_VERSION}
    && ${GZIP} -d ${PACKAGE_NAME}-${PACKAGE_VERSION}.tar.gz
    && ${TAR} -xf ${PACKAGE_NAME}-${PACKAGE_VERSION}.tar
    && cd ${PACKAGE_NAME}-${PACKAGE_VERSION}/
    && chmod u+w . && mkdir -p _build && mkdir -p _inst
    && chmod u+rwx _build _inst && chmod a-w .
    && cat ${CMAKE_BINARY_DIR}/CMakeCache.txt | ${SED} -n '/CMAKE_CACHEFILE_DIR:INTERNAL./{n\;x\;d\;}\;x\;1d\;\$G\;p'
       || ${SED} -n '/CMAKE_HOME_DIRECTORY:INTERNAL./{n\;x\;d\;}\;x\;1d\;\$G\;p'
       || ${SED} -n '/CMAKE_CXX_COMPILER:FILEPATH./{n\;x\;d\;}\;x\;1d\;\$G\;p'
       || ${SED} -n '/CMAKE_CXX_FLAGS:STRING./{n\;x\;d\;}\;x\;1d\;\$G\;p'
       || ${SED} -n '/CMAKE_CXX_FLAGS_DEBUG:STRING./{n\;x\;d\;}\;x\;1d\;\$G\;p'
       || ${SED} -n '/CMAKE_CXX_FLAGS_MINSIZEREL:STRING./{n\;x\;d\;}\;x\;1d\;\$G\;p'
       || ${SED} -n '/CMAKE_CXX_FLAGS_RELEASE:STRING./{n\;x\;d\;}\;x\;1d\;\$G\;p'
       || ${SED} -n '/CMAKE_CXX_FLAGS_RELWITHDEBINFO:STRING./{n\;x\;d\;}\;x\;1d\;\$G\;p'
       || ${SED} -n '/CMAKE_CXX_COMPILER-ADVANCED:INTERNAL./{n\;x\;d\;}\;x\;1d\;\$G\;p'
       || ${SED} -n '/CMAKE_CXX_COMPILER_WORKS:INTERNAL./{n\;x\;d\;}\;x\;1d\;\$G\;p'
       || ${SED} -n '/CMAKE_CXX_FLAGS-ADVANCED:INTERNAL./{n\;x\;d\;}\;x\;1d\;\$G\;p'
       || ${SED} -n '/CMAKE_CXX_FLAGS_DEBUG-ADVANCED:INTERNAL./{n\;x\;d\;}\;x\;1d\;\$G\;p'
       || ${SED} -n '/CMAKE_CXX_FLAGS_MINSIZEREL-ADVANCED:INTERNAL./{n\;x\;d\;}\;x\;1d\;\$G\;p'
       || ${SED} -n '/CMAKE_CXX_FLAGS_RELEASE-ADVANCED:INTERNAL./{n\;x\;d\;}\;x\;1d\;\$G\;p'
       || ${SED} -n '/CMAKE_CXX_FLAGS_RELWITHDEBINFO-ADVANCED:INTERNAL./{n\;x\;d\;}\;x\;1d\;\$G\;p'
       || ${SED} -n '/CMAKE_DETERMINE_CXX_ABI_COMPILED:INTERNAL./{n\;x\;d\;}\;x\;1d\;\$G\;p'
       || ${SED} -n '/CMAKE_C_COMPILER:FILEPATH./{n\;x\;d\;}\;x\;1d\;\$G\;p'
       || ${SED} -n '/CMAKE_C_FLAGS:STRING./{n\;x\;d\;}\;x\;1d\;\$G\;p'
       || ${SED} -n '/CMAKE_C_FLAGS_DEBUG:STRING./{n\;x\;d\;}\;x\;1d\;\$G\;p'
       || ${SED} -n '/CMAKE_C_FLAGS_MINSIZEREL:STRING./{n\;x\;d\;}\;x\;1d\;\$G\;p'
       || ${SED} -n '/CMAKE_C_FLAGS_RELEASE:STRING./{n\;x\;d\;}\;x\;1d\;\$G\;p'
       || ${SED} -n '/CMAKE_C_FLAGS_RELWITHDEBINFO:STRING./{n\;x\;d\;}\;x\;1d\;\$G\;p'
       || ${SED} -n '/CMAKE_C_COMPILER-ADVANCED:INTERNAL./{n\;x\;d\;}\;x\;1d\;\$G\;p'
       || ${SED} -n '/CMAKE_C_FLAGS-ADVANCED:INTERNAL./{n\;x\;d\;}\;x\;1d\;\$G\;p'
       || ${SED} -n '/CMAKE_C_FLAGS_DEBUG-ADVANCED:INTERNAL./{n\;x\;d\;}\;x\;1d\;\$G\;p'
       || ${SED} -n '/CMAKE_C_FLAGS_MINSIZEREL-ADVANCED:INTERNAL./{n\;x\;d\;}\;x\;1d\;\$G\;p'
       || ${SED} -n '/CMAKE_C_FLAGS_RELEASE-ADVANCED:INTERNAL./{n\;x\;d\;}\;x\;1d\;\$G\;p'
       || ${SED} -n '/CMAKE_C_FLAGS_RELWITHDEBINFO-ADVANCED:INTERNAL./{n\;x\;d\;}\;x\;1d\;\$G\;p'
       || ${SED} -n '/CMAKE_DETERMINE_C_ABI_COMPILED:INTERNAL./{n\;x\;d\;}\;x\;1d\;\$G\;p'
       || ${SED} -n '/EXECUTABLE_OUTPUT_PATH:PATH./{n\;x\;d\;}\;x\;1d\;\$G\;p'
       || ${SED} -n '/LIBRARY_OUTPUT_PATH:PATH./{n\;x\;d\;}\;x\;1d\;\$G\;p'
       || ${SED} -n '/produce slightly less optimized./{n\;x\;d\;}\;x\;1d\;\$G\;p'
       || cat -s > _build/CMakeCache.txt
    && cd _build
    && cmake -DCMAKE_INSTALL_PREFIX=${INSTDIR} .. || cmake ..
        || (echo "ERROR: the cmake configuration failed." && false)
    && make -j4
        || (echo "ERROR: the compilation failed." && false)
    && make test
        || (echo "ERROR: the test suite failed." && false)
    && make install
        || (echo "ERROR: the install target failed." && false)
    && make uninstall
        || (echo "ERROR: the uninstall target failed." && false)
    && test x`find ${INSTDIR} -type f | wc -l` = x0
        || (echo "ERROR: the uninstall target does not work." && false)
    && make clean
        || (echo "ERROR: the clean target failed." && false)
    && cd ${CMAKE_BINARY_DIR}/${PACKAGE_NAME}-${PACKAGE_VERSION}
    && chmod u+w . _build _inst && rm -rf _build _inst
    && find . -type d -print0 | xargs -0 chmod u+w
    && echo "=============================================================="
    && echo "${PACKAGE_NAME}-${PACKAGE_VERSION}"
            "is ready for distribution."
    && echo "=============================================================="
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Checking generated tarball..."
    )
  ADD_DEPENDENCIES(distcheck dist)
ENDMACRO()

