# Accurate mapping of tRNA reads
Refer to the original workflow: https://github.com/AnneHoffmann/tRNA-read-mapping

The mapping pipeline in this repository is identical, except for the software versions used. 
In addition to the multimapper protocols in the original code, we extend the pipeline for isoacceptor quantification. Once reads are mapped, all reads that map unambigously to a certain isoacceptor are considered for quantification (even though they might be ambiguously mapping to several isodecoders that share the same anticodon loop). This protocol is added as the variable:
multimapperHandling="quant"

For questions and help, please contact: xavier.hernandez@crg.eu
