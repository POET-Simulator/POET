simple_named_vec <- function(input) {
  input["Ca"] <- 0

  return(input)
}

extended_named_vec <- function(input, additional) {
  return(input["H"] + additional["H"])
}

bool_named_vec <- function(input) {
  return(input["H"] == 0)
}

simple_field <- function(field) {
  field <- as.data.frame(field, check.names = FALSE)
  field$Na <- 0
  return(field)
}

extended_field <- function(field, additional) {
  field <- as.data.frame(field, check.names = FALSE)
  additional <- as.data.frame(additional, check.names = FALSE)
  return(field + additional)
}
