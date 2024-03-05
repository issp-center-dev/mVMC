add_subdirectory(en/source)
add_subdirectory(ja/source)
add_subdirectory(tutorial/ja/source)

add_custom_target(doc DEPENDS doc-ja doc-en)
add_custom_target(tutorial DEPENDS tutorial-ja)
