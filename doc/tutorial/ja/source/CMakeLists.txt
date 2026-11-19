find_program(SPHINX_EXECUTABLE NAMES sphinx-build
  HINTS $ENV{SPHINX_DIR}
  PATH_SUFFIXES bin
  )

if(NOT SPHINX_EXECUTABLE)
  message(FATAL_ERROR "SPHINX_EXECUTABLE (sphinx-build) is not found!")
endif()

set(SPHINX_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR}/_build)
set(SPHINX_CACHE_DIR ${CMAKE_CURRENT_BINARY_DIR}/_doctrees)
set(SPHINX_HTML_DIR ${CMAKE_CURRENT_BINARY_DIR}/html)
set(SPHINX_PDF_DIR ${CMAKE_CURRENT_BINARY_DIR}/pdf)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/conf.py ${CMAKE_CURRENT_BINARY_DIR})

add_custom_target(tutorial-ja-html ALL
  COMMAND
    ${SPHINX_EXECUTABLE}
    ${SPHINX_OPTS}
    -b html
    -d ${SPHINX_CACHE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${SPHINX_HTML_DIR}
  )

add_custom_target(tutorial-ja-pdf ALL
  COMMAND
    ${SPHINX_EXECUTABLE}
    ${SPHINX_OPTS}
    -b latex
    -d ${SPHINX_CACHE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${SPHINX_PDF_DIR}
  # Without these options LaTeX stops at the first error and waits for
  # keyboard input. Not batchmode: the error must still reach the console.
  COMMAND
    ${CMAKE_COMMAND} -E env
    "LATEXOPTS=-interaction=nonstopmode -halt-on-error"
    make -C ${SPHINX_PDF_DIR}
  VERBATIM
  )

add_custom_target(tutorial-ja DEPENDS tutorial-ja-html tutorial-ja-pdf)

install(DIRECTORY ${SPHINX_HTML_DIR}
  DESTINATION doc/tutorial/ja
  PATTERN "_sources" EXCLUDE
  PATTERN ".buildinfo" EXCLUDE
  PATTERN "objects.inv" EXCLUDE
  )
install(DIRECTORY ${SPHINX_PDF_DIR}
  DESTINATION doc/tutorial/ja
  FILES_MATCHING PATTERN "*.pdf"
  )
