Add tsv files to the `mouse_mat-bl6` and `mouse_pat-cast` directories with the following format for step2.

| gene_id            | gene_name | chromosome | count |
| ------------------ | --------- | ---------- | ----- |
| ENSMUSG00000059327 | Eda       | X          | 10    |
| ENSMUSG00000009670 | Tex11     | X          | 14    |

Group the data for each mouse, one per directory. Please follow the following naming convention for your files:

```
(paternal or maternal)-(bl6 or cast)-(mouse_name).tsv
```

Example: `paternal-cast-mouse_1.tsv`