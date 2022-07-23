" For when editing one of test/test*.py
" Run one
map \f :w<C-m>:!clear; python -m pytest --capture=tee-sys %<C-m>
" Run all
map \t :w<C-m>:!clear; python -m pytest<C-m>
