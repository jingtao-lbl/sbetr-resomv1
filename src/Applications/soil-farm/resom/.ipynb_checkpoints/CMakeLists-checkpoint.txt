add_subdirectory(resomPara)
add_subdirectory(resom1layer)
add_subdirectory(resomNlayer)

set(BETR_LIBRARIES resomPara resom1layer resomNlayer;${BETR_LIBRARIES} PARENT_SCOPE)
