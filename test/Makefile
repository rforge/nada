PREFIX='/home/rclee/local/src/R-2.14.1'
R=${PREFIX}/bin/R
PACKAGE_DIR=../pkg
VERSION=`grep 'Version' ${PACKAGE_DIR}/DESCRIPTION \
         |cut -d':' -f2|tr -d '[:blank:]'`
PACKAGE=NADA
PACKAGE_FILE=${PACKAGE}_${VERSION}.tar.gz

default:
	echo "Usage make [build|remove|install|check|clean]"

build:
	${R} CMD build ${PACKAGE_DIR}

remove:
	${R} CMD REMOVE ${PACKAGE}

install: build
	${R} CMD INSTALL ${PACKAGE_FILE}

check: build
	${R} CMD check ${PACKAGE_FILE}

clean:
	rm -f ${PACKAGE_FILE}
	rm -rf ${PACKAGE}.Rcheck

