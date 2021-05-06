functions{
  real matrix_normal_lpdf(matrix X, matrix M, matrix U, matrix V) {
      int n = rows(X);
      int p = cols(X);
      // trace(D * B' * A * B)
      return -0.5 * trace_gen_quad_form(inverse(V), inverse(U), X - M) 
        - n*p/2.0 * log(2 * pi()) 
        - n/2.0 * log(determinant(V)) 
        - p/2.0 * log(determinant(U));
    }
}
