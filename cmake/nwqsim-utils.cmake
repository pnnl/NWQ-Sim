function(nwqsim_is_valid __variable __out)
  set(${__out} FALSE PARENT_SCOPE)
  if(DEFINED ${__variable} AND (NOT "${${__variable}}" STREQUAL ""))
      set(${__out} TRUE PARENT_SCOPE)
  endif()
endfunction()

#
# Sets an option's value if the user doesn't supply one.
#
function(nwqsim_option name value)
    nwqsim_is_valid(${name} was_set)
    if(was_set)
        message(STATUS "Value of ${name} was set by user to : ${${name}}")
    else()
        set(${name} ${value} PARENT_SCOPE)
        message(STATUS "Setting value of ${name} to default : ${value}")
    endif()
endfunction()