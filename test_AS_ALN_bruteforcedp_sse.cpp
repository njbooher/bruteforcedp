#include <cstdio>
#include "AS_global.H"
#include "AS_ALN_bruteforcedp_sse.H"

int main() {

  char *seq_a = "ACACACTAACACACTA";
  char *seq_b = "AGCACACAAGCACACA";
  
  alignLinker_s al;
  //dpCell M[AS_READ_MAX_NORMAL_LEN + 1][AS_READ_MAX_NORMAL_LEN + 1];

  dpCell *M = (dpCell *) safe_calloc((AS_READ_MAX_NORMAL_LEN + 1) * (AS_READ_MAX_NORMAL_LEN + 1), sizeof(dpCell));
  char *h_alignA = (char*) safe_calloc(2 * AS_READ_MAX_NORMAL_LEN + 2, sizeof(char));
  char *h_alignB = (char*) safe_calloc(2 * AS_READ_MAX_NORMAL_LEN + 2, sizeof(char));
  for (int i = 0; i < 50000; i++) {
    alignLinker_s aa;
    
    alignLinker(h_alignA,
                h_alignB,
                seq_a,
                seq_b,
                (dpCell (*)[AS_READ_MAX_NORMAL_LEN + 1]) M,
                &aa);
    al = aa;
  }
  
  printf("matches: %d\n", al.matches);
  printf("alignLen: %d\n", al.alignLen);
  printf("begI: %d\n", al.begI);
  printf("begJ: %d\n", al.begJ);
  printf("endI: %d\n", al.endI);
  printf("endJ: %d\n", al.endJ);
  printf("lenA: %d\n", al.lenA);
  printf("lenB: %d\n", al.lenB);
  printf("pIdentity: %f\n", al.pIdentity);
  printf("pCoverageA: %f\n", al.pCoverageA);
  printf("pCoverageB: %f\n", al.pCoverageB);

  safe_free2(h_alignB);
  safe_free2(h_alignA);
  safe_free2(M);

  return 0;
}
