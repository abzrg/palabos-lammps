# Function to check if an environment variable is defined
function(set_env_var_or_default var_name default_value)
    if(DEFINED ENV{${var_name}})
        set(${var_name} "$ENV{${var_name}}" PARENT_SCOPE)
    else()
        set(${var_name} "${default_value}" PARENT_SCOPE)
    endif()
    message(STATUS "${var_name} is set to ${${var_name}}")
endfunction()
