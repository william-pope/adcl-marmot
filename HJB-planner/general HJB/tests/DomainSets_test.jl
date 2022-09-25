using DomainSets: ×

# ?: is there any way to define a rectangle without using the cross symbol?
# my_rect = (0..2.1) × (0..3.4)

my_rect = Rectangle((0..2.8), (0..3.5))

@show typeof(my_rect)
@show eltype(my_rect)

@show my_rect.a
@show my_rect.b
@show my_rect.b[1]