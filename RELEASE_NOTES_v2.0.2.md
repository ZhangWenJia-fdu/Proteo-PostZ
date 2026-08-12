# ProteoPostZ Formal V2.0.2 Release Notes

ProteoPostZ Formal V2.0.2 updates PEAKS Online protein-result input while preserving the V2.0.1 analysis framework and other importers.

## PEAKS input

- Upload one PEAKS protein result file. The app detects DB versus LFQ from the internal column schema, not from the filename.
- DB schema uses `Area <sample>` for quantitative values and sample-specific `Coverage(%) <sample> > 0` for identification evidence.
- LFQ schema uses `<sample> Area` for quantitative values and excludes `Group N Area`, profile/ratio fields, significance, and other non-sample fields.
- If DB and LFQ-style fields are mixed, DB schema takes precedence.
- Both subtypes retain the first `Top == TRUE` row per `Protein Group`; if no Top row exists, the first source row in the group is retained. No within-group Area aggregation is performed.

## Count semantics

- PEAKS DB: `Identified_Protein_Count` uses sample-specific Coverage > 0; `Quantified_Protein_Count` uses non-missing Area.
- PEAKS LFQ: `Identified_Protein_Count` is unavailable (`NA`); `Quantified_Protein_Count` uses non-missing LFQ Area.
- Existing FragPipe/MSFragger, MaxQuant, DIA-NN, Spectronaut, and standard-matrix definitions remain unchanged.

## Runtime

- The Windows package includes the portable R runtime, portable R libraries, app source, annotations, launcher, scripts, and bilingual instructions.
- Normal operation is offline and opens the local Shiny interface at `http://127.0.0.1:3840/`.
- Local test files, test records, logs, generated outputs, and source-repository history are not included in the package.
