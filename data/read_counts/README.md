Add tsv files to the `mouse_mat` and `mouse_pat` directories with the following format for step0.

| gene_id            | gene_name | chromosome | count |
| ------------------ | --------- | ---------- | ----- |
| ENSMUSG00000059327 | Eda       | X          | 10    |
| ENSMUSG00000009670 | Tex11     | X          | 14    |

Group the data for each mouse, one per directory. Please follow the following naming convention for your files:

```bash
(mouse_id)-(gender of mouse)-(chromosome parental origin: pat or mat)-(strain chromosome: bl6 or cast).tsv
```

Example: `mouse_1-female-pat-cast.tsv`

It is important that the files are named this way, because the information is used during the pipeline to distinguish between samples.