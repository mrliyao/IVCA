* This solution follows the below work and project:

  > Menon, V. V., Feldmann, C., Amirpour, H., Ghanbari, M., & Timmerer, C. (2022, June). VCA: video complexity analyzer. In *Proceedings of the 13th ACM multimedia systems conference* (pp. 259-264).
  >
  > VCA: Video Complexity Analyzer, https://github.com/cd-athena/VCA

* Usage
  1. Test dataset directory needs to be set in `submisssion.sh`
  2. The resolution and chroma sub-sampling format of the test sequence need to be set in `submission.sh`, per sequence.
  3. Default target platform is Linux with amd64 architecture, using the execution  `vca` and `sum`, Windows with amd64 architecture is also supported with the execution  `vca.exe` and `sum.exe`
  4. For building the the execution `vca` and `sum` respectively, please refer to `./source/VCA/` and `./source/sum.cpp`, the former provides a build guide while the latter can be build with `g++ ./sum.cpp -o sum`

* Contact
  * Yao Li, mrliyao@mail.ustc.edu.cn
  * Junqi Liao, liaojq@mail.ustc.edu.cn
  * Zhuoyuan Li, zhuoyuanli@mail.ustc.edu.cn
  * Advisor: Prof. Li Li, lil1@ustc.edu.cn and Prof. Dong Liu, dongeliu@ustc.edu.cn

