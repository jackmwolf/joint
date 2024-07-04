# When performing super learning based model averaging, joint requires
# the Rsolnp package for metalearning. This package is only suggested by
# dependency sl3 and an error will be generated if Rsolnp is not installed.

# This workaround avoids a NOTE about Rsolnp being imported but not used.

ignore_unused_imports <- function() {
  Rsolnp::solnp
}
