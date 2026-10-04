## Test environments
* local macOS 15.7.9 M1 install, R 4.5.2
* devtools::check()
* devtools::release_checks()
* devtools::check_mac_release()
* rhub::rhub_check()
* https://win-builder.r-project.org/upload.aspx

## R CMD check results
There were no ERRORs, or WARNINGs. There was 1 NOTE. Win-builder flags my email address 
again? It also mentions an invalid URL, but when I copy paste it, it does work. 

## Downstream dependencies
* revdepcheck::revdep_check(num_workers = 4)

No problems found, still no other packages depends on this package
