# lasso_extensions

## English description

This repository contains the project for the course Statistical Methods for High Dimensional Data. We needed to replicate the simulation study presented in Laurent et al. [2009](https://dl.acm.org/doi/pdf/10.1145/1553374.1553431?casa_token=q0cqNqvBgU8AAAAA:3UqYS8uH6Z4OEYBHWctGsPv7lS0S748lov5D4B-u4SmGt9WcF9lFwTdVwRYfdDmmHZvdSpbxS08), as well as the simulation study and the real data application presented in Lim & Hastie [2014](https://hastie.su.domains/Papers/glinternet_jcgs.pdf).

For the first article we needed to show on simulated data how overlapping group lasso greatly outperfors group lasso when the covariate groups are actually overlapping. For the second one, first we had to asses how the proposed approach, a sparse linear model with interactions estimed via overlapping group lasso penalty, was able to retrieve real intercation terms using simulated data. Secondly, we used real-world data to compare the proposed model to a similar technique and to a boosted tree with first-order interactions.

## Italian description

Questa repository contiene il progetto per il corso di Metodi Statistici per Dati ad Elevata. La richiesta era quella di replicare lo studio di simulazione presentato in Laurent et al. [2009](https://dl.acm.org/doi/pdf/10.1145/1553374.1553431?casa_token=q0cqNqvBgU8AAAAA:3UqYS8uH6Z4OEYBHWctGsPv7lS0S748lov5D4B-u4SmGt9WcF9lFwTdVwRYfdDmmHZvdSpbxS08), come anche lo studio di simulazione e l'applicazione su dati reali presenti in Lim & Hastie [2014](https://hastie.su.domains/Papers/glinternet_jcgs.pdf).

Per il primo articolo dovevamo mostrare su dati simulati come l'overlap group lasso superi notevolmente il group lasso quando i gruppi di covariate si sovrappongono. Per il secondo, in primo luogo abbiamo dovuto valutare come l'approccio proposto, un modello lineare sparso con interazioni stimate tramite penalità di tipo overlap group lasso, fosse in grado di identificare i veri termini di interazione, utilizzando dati simulati. In secondo luogo, abbiamo utilizzato i dati del mondo reale per confrontare il modello proposto con una tecnica simile e con un boosting con interazioni del primo ordine.
