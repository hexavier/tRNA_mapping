# Accurate mapping of tRNA reads
Refer to the original workflow for a detailed description of the pipeline and a demo dataset: 
https://github.com/AnneHoffmann/tRNA-read-mapping

Updates with regard to the original workflow:
- The mapping pipeline in this repository is identical, except for the software versions used. All versions and details are commented within the "best_practice_workflow.sh".
- In addition to the multimapper protocols in the original code, we extend the pipeline for isoacceptor quantification. Once reads are mapped, all reads that map unambigously to a certain isoacceptor are considered for quantification (even though they might be ambiguously mapping to several isodecoders that share the same anticodon loop). This protocol is added as the variable:
multimapperHandling="quant"
- Scripts have been adapted to run on a Sun Grid Engine (SGE) cluster using the qsub command.

For questions and help, please contact: xavier.hernandez@crg.eu
