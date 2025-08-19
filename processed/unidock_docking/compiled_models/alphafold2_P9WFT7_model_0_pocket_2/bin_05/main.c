
#include "header.h"


static const int32_t num_class[] = {  2, };

int32_t get_num_target(void) {
  return N_TARGET;
}
void get_num_class(int32_t* out) {
  for (int i = 0; i < N_TARGET; ++i) {
    out[i] = num_class[i];
  }
}
int32_t get_num_feature(void) {
  return 1940;
}
const char* get_threshold_type(void) {
  return "float64";
}
const char* get_leaf_output_type(void) {
  return "float64";
}

void predict(union Entry* data, int pred_margin, double* result) {
  unsigned int tmp;
  predict_unit0(data, result);
  predict_unit1(data, result);
  predict_unit2(data, result);
  predict_unit3(data, result);
  predict_unit4(data, result);
  predict_unit5(data, result);
  predict_unit6(data, result);
  predict_unit7(data, result);
  predict_unit8(data, result);
  predict_unit9(data, result);
  predict_unit10(data, result);
  predict_unit11(data, result);
  predict_unit12(data, result);
  predict_unit13(data, result);
  predict_unit14(data, result);
  predict_unit15(data, result);
  
  // Average tree outputs
  result[0] /= 501;
  result[1] /= 501;
  
  // Apply base_scores
  result[0] += 0;
  result[1] += 0;
  
  // Apply postprocessor
  if (!pred_margin) { postprocess(result); }
}

void postprocess(double* result) {
  // Do nothing
}

