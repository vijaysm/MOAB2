This debian folder tries to support various Debian/Ubuntu releases.

Build
-----
Simply do
```
make
```
and the respective `control` and `rules` will be built.

Unfortunately, pre-Wily Ubuntu releases have a Metis bug that makes it
impossible to be used with MOAB. For this reason, MOAB build without Metis
support by default.

If you want to build with Metis or any other additional dependency that cannot
be included by default, use the environment variariables `ADDITIONAL_DEPS` and
`ADDITIONAL_ENABLES`. To build MOAB with support for Metis and ParMetis, for
example, do
```
ADDITIONAL_DEPS=", libmetis-dev, libparmetis-dev" \
ADDITIONAL_ENABLES="-DENABLE_METIS:BOOL=ON -DENABLE_PARMETIS:BOOL=ON" \
make
```
