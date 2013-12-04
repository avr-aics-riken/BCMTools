105,106c105,106
<   comm.Allreduce(MPI_IN_PLACE, &levelMin, 1, MPI::INTEGER, MPI::MIN);
<   comm.Allreduce(MPI_IN_PLACE, &levelMax, 1, MPI::INTEGER, MPI::MAX);
---
>   comm.Allreduce(MPI_IN_PLACE, &levelMin, 1, MPI::INT, MPI::MIN);
>   comm.Allreduce(MPI_IN_PLACE, &levelMax, 1, MPI::INT, MPI::MAX);
127c127
<   comm.Reduce(&numBlock, &numBlockSum, 1, MPI::INTEGER, MPI::SUM, 0);
---
>   comm.Reduce(&numBlock, &numBlockSum, 1, MPI::INT, MPI::SUM, 0);
129c129
<     comm.Reduce(MPI_IN_PLACE, nBlock, nLevel, MPI::INTEGER, MPI::SUM, 0);
---
>     comm.Reduce(MPI_IN_PLACE, nBlock, nLevel, MPI::INT, MPI::SUM, 0);
131c131
<     comm.Reduce(nBlock, nBlock, nLevel, MPI::INTEGER, MPI::SUM, 0);
---
>     comm.Reduce(nBlock, nBlock, nLevel, MPI::INT, MPI::SUM, 0);
145,146c145,146
<   comm.Reduce(&numBlock, &numBlockMin, 1, MPI::INTEGER, MPI::MIN, 0);
<   comm.Reduce(&numBlock, &numBlockMax, 1, MPI::INTEGER, MPI::MAX, 0);
---
>   comm.Reduce(&numBlock, &numBlockMin, 1, MPI::INT, MPI::MIN, 0);
>   comm.Reduce(&numBlock, &numBlockMax, 1, MPI::INT, MPI::MAX, 0);
149c149
<     comm.Reduce(MPI_IN_PLACE, &numBlock2Sum, 1, MPI::INTEGER, MPI::SUM, 0);
---
>     comm.Reduce(MPI_IN_PLACE, &numBlock2Sum, 1, MPI::INT, MPI::SUM, 0);
151c151
<     comm.Reduce(&numBlock2Sum, &numBlock2Sum, 1, MPI::INTEGER, MPI::SUM, 0);
---
>     comm.Reduce(&numBlock2Sum, &numBlock2Sum, 1, MPI::INT, MPI::SUM, 0);
211,212c211,212
<   comm.Reduce(nFaceInter, nFaceInterSum, 3, MPI::INTEGER, MPI::SUM, 0);
<   comm.Reduce(nFaceIntra, nFaceIntraSum, 3, MPI::INTEGER, MPI::SUM, 0);
---
>   comm.Reduce(nFaceInter, nFaceInterSum, 3, MPI::INT, MPI::SUM, 0);
>   comm.Reduce(nFaceIntra, nFaceIntraSum, 3, MPI::INT, MPI::SUM, 0);
235,236c235,236
<   comm.Reduce(&nFaceInterTotal, &nFaceInterTotalMin, 1, MPI::INTEGER, MPI::MIN, 0);
<   comm.Reduce(&nFaceInterTotal, &nFaceInterTotalMax, 1, MPI::INTEGER, MPI::MAX, 0);
---
>   comm.Reduce(&nFaceInterTotal, &nFaceInterTotalMin, 1, MPI::INT, MPI::MIN, 0);
>   comm.Reduce(&nFaceInterTotal, &nFaceInterTotalMax, 1, MPI::INT, MPI::MAX, 0);
240c240
<     comm.Reduce(MPI_IN_PLACE, &nFaceInterTotal2Sum, 1, MPI::INTEGER, MPI::SUM, 0);
---
>     comm.Reduce(MPI_IN_PLACE, &nFaceInterTotal2Sum, 1, MPI::INT, MPI::SUM, 0);
242c242
<     comm.Reduce(&nFaceInterTotal2Sum, &nFaceInterTotal2Sum, 1, MPI::INTEGER, MPI::SUM, 0);
---
>     comm.Reduce(&nFaceInterTotal2Sum, &nFaceInterTotal2Sum, 1, MPI::INT, MPI::SUM, 0);
