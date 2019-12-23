# sums_of_three_cubes

A non-official python implement of paper ["CRACKING THE PROBLEM WITH 33"](https://people.maths.bris.ac.uk/~maarb/papers/cubesv1.pdf) which present a computational solution to sums of three cubes problems for 33.

Now this project is just a toy python implement to the paper above. However, it accomplish the gist. Improvement may or may not plan to be done.

现在的实现还远不够完美，基本还处于玩具阶段，但基本实现了大体解题思路，还有一些没有实现的部分如下：

1. [Montgomery‘s trick](https://en.wikipedia.org/wiki/Montgomery_modular_multiplication),这个我怀疑python的pow函数就有内置这个算法;
2. M'; 由于M'的选用有很大的任意性，论文并没有给出M'的选用方法，所以这部分暂时没有去实现;
3. 算法根据d与$\sqrt {\alpha B}$的关系分成2部分;这样的话大于$\sqrt {\alpha B}$的素数的三次方表格就不需要计算了;
4. 多进程计算;
5. C实现;这个不打算做，主要是python写的惯;
6. 多进程;程序结构本质是多进程的，有空弄下;
7. 日志;现在居然一点日志都没有
8. 表格的组成;现在很多都弄成表格，速度快了，但内存占用太大，计算不了太大的值;因此可能考虑放开一些表格，实时但重复计算，做时间-空间交易;

大家可以用http://www.asahi-net.or.jp/~KC2H-MSM/mathland/math04/matb0100.htm找一些测试案例，但请注意，论文的方法只对$k\equiv3(mod\;9)$成立;

