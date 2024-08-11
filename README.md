# chem-env

`chem-env` provides Python tools for chemistry that can easily be deployed as modal apps.

## Usage

### Chemistry tools

### Modal app

Deploy the app by running

```
chemenv deploy
```

Once the app is deployed, you can run functions like

```
app_name = 'chemenv'
fxn = modal.Function.lookup(app_name, "get_tanimoto_similarity")
result = await fxn.remote.aio("CCO", "CC")
```
