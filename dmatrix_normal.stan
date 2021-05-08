functions{
  real dmatrix_normal(matrix Y, matrix M, matrix U, matrix V) {
      int n = rows(Y);
      int p = cols(Y);
      // trace(D * B' * A * B)
      return -0.5 * trace_gen_quad_form(inverse(V), inverse(U), Y - M) 
        - n*p/2.0 * log(2 * pi()) 
        - n/2.0 * log(determinant(V)) 
        - p/2.0 * log(determinant(U));
    }
}
