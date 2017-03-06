std::vector<std::vector<long double>> lagrange_nodes(std::vector<std::vector<long double>> u,std::vector<long double> q_weights) {
  for (unsigned int i=0;i<q_weights.size();i++) {
      for (unsigned int j=0;j<q_weights.size();j++) {
          u[i][j]*=q_weights[i]*q_weights[j];
        }
    }
  return u;
}
